"""
    OptResults

Outputs of one design evaluation ([`solve_structure`](@ref)), carrying
everything downstream objectives and constraints consume — all plain data,
so gradients of any scalar of these fields flow back to the design vector.

# Fields
- `U::Vector`: full-space displacements, 6 slots per node in the order
  `ux, uy, uz, θx, θy, θz` — node `i`'s translations are `U[6(i−1) .+ (1:3)]`
  (e.g. all z-displacements are `U[3:6:end]`)
- `X::Matrix`: evaluated node positions (3 × n_nodes)
- `A::Vector`: evaluated per-element areas
- `L::Vector`: evaluated element lengths
- `sections::AbstractVector`: the evaluated per-element sections (feed to
  further pure computations if needed)
- `EA::AbstractVector`: evaluated axial rigidities (plain vector — cheap
  under reverse AD)
- `F::Vector`: the full-space load vector the design was solved for —
  `P − Pf(x)` with design-dependent element loads (`p.loads = :design`),
  else the base-design `p.F`
- `ends`: the evaluated per-element `EndConditions` (joint variables), or
  `nothing` when the problem has none
"""
struct OptResults{TU,TX,TA,TL,TS,TE,TF,TN}
    U::TU
    X::TX
    A::TA
    L::TL
    sections::TS
    EA::TE
    F::TF
    ends::TN
end

"""
    _design_state(x, p::OptParams) -> (X, A, sections, EAvec, ends)

The shared pure core of a design evaluation: scatter the design vector into
node positions (`X0 + Sx·x`, additive), per-element areas (`Sa·x` where
masked, `A0` elsewhere), rebuilt sections, axial rigidities `E·A`, and
joint-substituted end conditions. Everything downstream —
[`solve_structure`](@ref) and [`updatemodel`](@ref) — is built on this one
function, so the optimization path and the materialized model can never
disagree.
"""
function _design_state(x::AbstractVector, p::OptParams)
    X = p.X0 + reshape(p.Sx * x, 3, :)
    # masked scatter as arithmetic BROADCASTS, not a comprehension (per-
    # element getindex pullbacks cost ~35×) and not `ifelse.` (its mixed
    # Dual/Float64 branches give a union eltype that poisons downstream
    # sparse algebra under ForwardDiff)
    A = p.A0 .* .!p.amask .+ (p.Sa * x) .* p.amask
    # no section-affecting variables ⇒ the sections are constants: reuse
    # them instead of rebuilding (and differentiating) identical structs.
    # Flexural/torsional properties evaluate only when some element has a
    # SectionVariable (props falls back to the constant base vector, and
    # for unmasked rows the two coincide) — the closure keeps the simple
    # two-branch mask shape, which Enzyme compiles where a three-branch
    # variant hits an internal error.
    secvar = any(p.pmask)
    props = secvar ? p.P0v .* .!p.prmask .+ (p.Sp * x) .* p.prmask : p.P0v
    secmask = p.amask .| p.pmask
    sections = if any(secmask)
        map(eachindex(secmask)) do i
            s = p.base_sections[i]::Section{Float64}
            secmask[i] ? Section(s.material, A[i], props[3i-2], props[3i-1], props[3i]) : s
        end
    else
        p.base_sections
    end
    if any(p.tmask)
        t = p.St * x
        sections = map(eachindex(p.tmask)) do i
            p.tmask[i] ? ShellSection(p.base_sections[i].material, t[i]) : sections[i]
        end
    end
    EAvec = p.Evec .* A
    ends = _design_ends(x, p)
    return X, A, sections, EAvec, ends
end

"""
    _design_ends(x, p::OptParams) -> Vector{EndConditions} | nothing

Per-element `EndConditions` with joint-variable rotational stiffnesses
(`ky = kz = factor·x[slot]`) substituted at the flagged ends. Returns
`nothing` when the problem has no joint variables — Asap's frame kernels
then use the model's cached ends, so joint-free problems pay zero overhead.
"""
function _design_ends(x::AbstractVector, p::OptParams)
    all(iszero, p.jslot1) && all(iszero, p.jslot2) && return nothing
    return map(eachindex(p.base_ends)) do i
        s1, s2 = p.jslot1[i], p.jslot2[i]
        s1 == 0 && s2 == 0 && return p.base_ends[i]
        base = p.base_ends[i]::EndConditions{Float64}
        e1 = s1 == 0 ? base.e1 :
             EndSprings(base.e1.kx, base.e1.kt, p.jfac1[i] * x[s1], p.jfac1[i] * x[s1])
        e2 = s2 == 0 ? base.e2 :
             EndSprings(base.e2.kx, base.e2.kt, p.jfac2[i] * x[s2], p.jfac2[i] * x[s2])
        EndConditions(e1, e2)
    end
end

"""
    _element_lengths(X, p::OptParams) -> Vector

Per-element Euclidean lengths from evaluated node positions `X` (3 × n):
column norms of the batched end-to-end difference `X · Cincᵀ`. One matmul +
one reduction — a handful of AD graph nodes total, where the per-element
formulation cost one pullback closure per element.
"""
function _element_lengths(X, p::OptParams)
    ΔX = X * transpose(p.Cinc)               # 3 × n_el
    return vec(sqrt.(sum(abs2, ΔX; dims=1)))
end

"""
    solve_structure(x, p::OptParams) -> OptResults

Evaluate a design vector: scatter positions and areas, build the load
vector (element-load fixed-end forces re-evaluated at the design when
`p.loads = :design`), run Asap's pure solve on the resulting `ModelState`,
and package the results. Fully differentiable —
`Zygote.gradient(x -> f(solve_structure(x, p)), x0)` needs nothing else, and
the gradient includes the load path.

`solve_truss` and `solve_frame` are aliases (the unified core needs no
distinction).
"""
function solve_structure(x::AbstractVector, p::OptParams)
    X, A, sections, EAvec, ends = _design_state(x, p)
    # eltype follows the design vector: Float64 normally, ForwardDiff Duals
    # under forward-mode AD (the Dual solve lives in AsapForwardDiffExt)
    state = ModelState{eltype(X)}(X, sections, EAvec, ends)
    cache = p.model.cache
    K = Asap.assemble_K(cache, state)
    free = cache.partition.free
    # two branches behind the Asap._pure_solve barrier: the constant-load
    # path keeps a concretely-typed solve call (Enzyme needs it); the
    # design-load path evaluates F(x) once and reuses it for compliance
    cm = cache.constraint
    if cm !== nothing   # rigid diaphragms: the reduced system CᵀKC·û = CᵀF (Asap ≥ 1.3)
        F = _design_F(p) ? p.P - Asap.assemble_Pf(cache, state) : p.F
        U = cache.free_embed * Asap._constrained_solve(cm, K, F[free], p.solver)
        return OptResults(U, X, A, _element_lengths(X, p), sections, EAvec, F, ends)
    end
    if _design_F(p)
        F = p.P - Asap.assemble_Pf(cache, state)
        U = cache.free_embed * Asap._pure_solve(K, F[free], p.solver)
        return OptResults(U, X, A, _element_lengths(X, p), sections, EAvec, F, ends)
    else
        U = cache.free_embed * Asap._pure_solve(K, p.F[free], p.solver)
        return OptResults(U, X, A, _element_lengths(X, p), sections, EAvec, p.F, ends)
    end
end

"does this problem re-evaluate element-load fixed-end forces per design?"
# FRAME element loads follow the design (gates the frame-specific corrections)
_design_loads(p::OptParams) = p.loads === :design && size(p.Eload, 2) > 0
# the load vector F follows the design: frame loads OR shell surface loads
# (self-weight ∝ ρt, pressure ∝ area/orientation — Asap ≥ 1.3)
_design_F(p::OptParams) = _design_loads(p) || (p.loads === :design && !isempty(p.model.cache.fef.sloads))

const solve_truss = solve_structure
const solve_frame = solve_structure

"""
    axial_force(res::OptResults, p::OptParams) -> Vector

Per-element axial force [force], tension-positive, matching Asap's stateful
`axial_force(res, el)` (slot 7 of the local end forces `f = k·u_local + q̃`:
the axial action at the END node): `EA/L` times the axial elongation, plus
the axial fixed-end force `q̃₇` of any element loads on the member (an
inclined loaded beam, self-weight on a column, …), evaluated at the design
when `p.loads = :design`. Pure and differentiable; exact for truss elements
and for frame elements with rigid axial connections.
"""
function axial_force(res::OptResults, p::OptParams)
    # batched: N = EA/L · (x̂ ⋅ Δu) = EA · (ΔX ⋅ ΔU) / L²  — two incidence
    # matmuls + broadcasts instead of one pullback closure per element
    ΔX = res.X * transpose(p.Cinc)                         # 3 × n_el
    ΔU = reshape(res.U, 6, :)[1:3, :] * transpose(p.Cinc)  # translational diffs
    N = vec(sum(ΔX .* ΔU; dims=1)) .* res.EA ./ (res.L .^ 2)
    size(p.Eload, 2) == 0 && return N
    _design_loads(p) || return N .+ p.qa0
    Q = Asap.fef_local(p.model.cache, ModelState{eltype(res.X)}(res.X, res.sections, res.EA, res.ends))
    return N .+ p.Eload * Q[7, :]
end

"""
    axial_stress(res::OptResults, p::OptParams) -> Vector

Per-element axial stress [force/length²]: axial force over evaluated area.
"""
axial_stress(res::OptResults, p::OptParams) = axial_force(res, p) ./ res.A

# sum∘broadcast, not LinearAlgebra.dot: dot lowers to BLAS, which Enzyme's
# forward mode cannot handle under runtime activity
"""
    compliance(res::OptResults, p::OptParams) -> Real

External work `Fᵀu` of the evaluated design — the canonical smooth
stiffness objective. `F` is the design's own load vector (`res.F`): with
design-dependent element loads it varies with the design, and the gradient
includes that term.
"""
compliance(res::OptResults, p::OptParams) = sum(res.F .* res.U)

"""
    GeometricProperties(x, p::OptParams)

Geometry-only evaluation of a design — positions, lengths, areas — with NO
structural solve. Use for objectives like material volume
(`dot(geo.A, geo.L)`) where paying for a solve would be waste.

# Fields
- `X::Matrix`, `L::Vector`, `A::Vector`
"""
struct GeometricProperties{TX,TL,TA}
    X::TX
    L::TL
    A::TA
end

function GeometricProperties(x::AbstractVector, p::OptParams)
    # geometry-only: skip section-struct construction entirely (its
    # constructor pullbacks would dominate this otherwise trivial path)
    X = p.X0 + reshape(p.Sx * x, 3, :)
    A = p.A0 .* .!p.amask .+ (p.Sa * x) .* p.amask
    return GeometricProperties(X, _element_lengths(X, p), A)
end

"""
    shell_volume(x, p::OptParams) -> T

Material volume of all `ShellSection` elements at the design `x`: Σ t·area,
with thickness from [`ThicknessVariable`](@ref)s (else the base thickness)
and area from the current node positions. Pure and differentiable — the
natural cost term for wall/slab thickness optimization.
"""
function shell_volume(x::AbstractVector, p::OptParams)
    isempty(p.sidx) && return zero(eltype(x))
    X = p.X0 + reshape(p.Sx * x, 3, :)
    t = any(p.tmask) ? p.T0 .* .!p.tmask .+ (p.St * x) .* p.tmask : p.T0
    return sum(map(eachindex(p.sidx)) do k
        c = p.snodes[:, k]
        xy = Asap.shell_frame(_col3(X, c[1]), _col3(X, c[2]), _col3(X, c[3]), _col3(X, c[4]))[2]
        t[p.sidx[k]] * Asap._q4_area(xy)
    end)
end

@inline _col3(X, i) = Asap.SVector(X[1, i], X[2, i], X[3, i])
