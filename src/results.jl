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
"""
struct OptResults{TU,TX,TA,TL,TS,TE}
    U::TU
    X::TX
    A::TA
    L::TL
    sections::TS
    EA::TE
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
    # masked scatter as a BROADCAST, not a comprehension — per-element
    # getindex pullbacks (one one-hot adjoint each) cost ~35× the
    # broadcasted ifelse under Zygote
    A = ifelse.(p.amask, p.Sa * x, p.A0)
    # no area variables ⇒ the sections are constants: reuse them instead of
    # rebuilding (and differentiating) an identical struct per element
    sections = if any(p.amask)
        map(eachindex(p.amask)) do i
            s = p.base_sections[i]::Section{Float64}
            p.amask[i] ? Section(s.material, A[i], s.Ix, s.Iy, s.J) : s
        end
    else
        p.base_sections
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

Evaluate a design vector: scatter positions and areas, run Asap's pure
solve on the resulting `ModelState`, and package the results. Fully
differentiable — `Zygote.gradient(x -> f(solve_structure(x, p)), x0)`
needs nothing else.

`solve_truss` and `solve_frame` are aliases (the unified core needs no
distinction).
"""
function solve_structure(x::AbstractVector, p::OptParams)
    X, A, sections, EAvec, ends = _design_state(x, p)
    state = ModelState{Float64}(X, sections, EAvec, ends)
    U = Asap.solve(p.model, state)
    return OptResults(U, X, A, _element_lengths(X, p), sections, EAvec)
end

const solve_truss = solve_structure
const solve_frame = solve_structure

"""
    axial_force(res::OptResults, p::OptParams) -> Vector

Per-element axial force [force], tension-positive: `EA/L` times the axial
elongation extracted from the global displacements. Pure and
differentiable; exact for truss elements and for the axial action of frame
elements with rigid axial connections.
"""
function axial_force(res::OptResults, p::OptParams)
    # batched: N = EA/L · (x̂ ⋅ Δu) = EA · (ΔX ⋅ ΔU) / L²  — two incidence
    # matmuls + broadcasts instead of one pullback closure per element
    ΔX = res.X * transpose(p.Cinc)                         # 3 × n_el
    ΔU = reshape(res.U, 6, :)[1:3, :] * transpose(p.Cinc)  # translational diffs
    return vec(sum(ΔX .* ΔU; dims=1)) .* res.EA ./ (res.L .^ 2)
end

"""
    axial_stress(res::OptResults, p::OptParams) -> Vector

Per-element axial stress [force/length²]: axial force over evaluated area.
"""
axial_stress(res::OptResults, p::OptParams) = axial_force(res, p) ./ res.A

"""
    compliance(res::OptResults, p::OptParams) -> Real

External work `Fᵀu` of the evaluated design — the canonical smooth
stiffness objective.
"""
compliance(res::OptResults, p::OptParams) = dot(p.F, res.U)

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
    A = ifelse.(p.amask, p.Sa * x, p.A0)
    return GeometricProperties(X, _element_lengths(X, p), A)
end
