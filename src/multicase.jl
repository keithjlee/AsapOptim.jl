"""
    LoadCases

Several load cases of one optimization problem, solved together against ONE
stiffness factorization per design (a multi-right-hand-side solve). Built
from an [`OptParams`](@ref) and the case tags of the model's loads:

    lc = LoadCases(p)                       # every case in the model, first-appearance order
    lc = LoadCases(p; cases = [:dead, :snow])

Everything here is constant per optimization; the design-dependent part of
each case's load vector (element-load fixed-end forces that follow member
geometry and sections, `p.loads = :design`) is re-evaluated per design by
[`solve_cases`](@ref).

A load combination is a linear combination of cases: solve the cases once
and combine `res.U * λ` (or build the combination as its own case by tagging
factored loads with one case name).

# Fields
- `names::Vector{Symbol}`: the case tags, one column each below
- `P::Matrix{Float64}`: full-space nodal loads per case (n_dof × n_cases)
- `F::Matrix{Float64}`: full-space `P − Pf` per case at the BASE design (the
  frozen load vectors)
- `static::Matrix{Float64}`: design-independent fixed-end forces per case
  (element loads on `VariableElement`s)
- `M::Matrix{Float64}`: 12·n_loads × n_cases — 1 where frame element load `j`
  (`cache.fef` order) belongs to the case, repeated over its 12 slots
- `E::Matrix{Float64}`: n_loads × n_cases case membership of each frame
  element load
- `qa0::Matrix{Float64}`: per-element axial fixed-end force `q̃₇` per case at
  the base design (frozen loads)
"""
struct LoadCases
    names::Vector{Symbol}
    P::Matrix{Float64}
    F::Matrix{Float64}
    static::Matrix{Float64}
    M::Matrix{Float64}
    E::Matrix{Float64}
    qa0::Matrix{Float64}
end

function LoadCases(p::OptParams; cases::AbstractVector{Symbol}=Asap.load_cases(p.model))
    isempty(cases) && throw(ArgumentError("no load cases"))
    allunique(cases) || throw(ArgumentError("load case names must be distinct"))
    model = p.model
    cache = model.cache
    present = Asap.load_cases(model)
    for c in cases
        c in present || throw(ArgumentError("no load in the model carries case :$c"))
    end
    t = cache.fef
    if p.loads === :design && !isempty(t.sloads)
        throw(ArgumentError("LoadCases with design-dependent SHELL surface loads is not supported yet " *
                            "(use loads = :frozen, or one OptParams per case)"))
    end
    nl = length(t.loads)
    nc = length(cases)
    E = zeros(nl, nc)
    for (j, l) in enumerate(t.loads), (k, c) in enumerate(cases)
        l.case === c && (E[j, k] = 1.0)
    end
    M = repeat(E; inner=(12, 1))
    ng = cache.partition.n_global
    nel = length(model.elements)
    P = zeros(ng, nc)
    F = zeros(ng, nc)
    S = zeros(ng, nc)
    qa0 = zeros(nel, nc)
    try
        for (k, c) in enumerate(cases)
            assemble_loads!(cache, model; case=c)
            P[:, k] = cache.P
            F[:, k] = cache.P .- cache.Pf
            S[:, k] = cache.fef.static
            for (e, el) in enumerate(model.elements)
                el isa FrameElement && (qa0[e, k] = cache.q_local[e][1][7])
            end
        end
    finally
        assemble_loads!(cache, model)          # restore the all-loads state the params rely on
    end
    return LoadCases(collect(cases), P, F, S, M, E, qa0)
end

Base.length(lc::LoadCases) = length(lc.names)

"""
    solve_cases(x, p::OptParams, lc::LoadCases) -> OptResults

Evaluate a design for every case in `lc` with one stiffness assembly and one
factorization: like [`solve_structure`](@ref), but `res.U` and `res.F` are
matrices (n_dof × n_cases, column `k` = case `lc.names[k]`). Positions,
sections and lengths are shared. Differentiable through the same paths
(ForwardDiff Duals, ChainRules engines — Asap's multi-RHS `solve_free`
rules).

Per-case scalars: [`compliance`](@ref)`(res, p, lc)`,
[`axial_force`](@ref)`(res, p, lc)`, [`axial_stress`](@ref)`(res, p, lc)`.

```julia
lc  = LoadCases(p)                       # e.g. :dead, :snow, :wind
res = solve_cases(x, p, lc)
N   = axial_force(res, p, lc)            # n_el × 3
Nc  = N * [1.0, 1.0, 0.0]                # D + S, by superposition
```
"""
function solve_cases(x::AbstractVector, p::OptParams, lc::LoadCases)
    X, A, sections, EAvec, ends = _design_state(x, p)
    state = ModelState{eltype(X)}(X, sections, EAvec, ends)
    cache = p.model.cache
    K = Asap.assemble_K(cache, state)
    free = cache.partition.free
    F = _cases_F(p, lc, state)
    cm = cache.constraint
    U = cm === nothing ? cache.free_embed * Asap._pure_solve(K, F[free, :], p.solver) :
        cache.free_embed * Asap._constrained_solve(cm, K, F[free, :], p.solver)
    return OptResults(U, X, A, _element_lengths(X, p), sections, EAvec, F, ends)
end

# per-case full-space load vectors: frozen, or P − (static + frame FEF of the
# case's loads at the design state)
function _cases_F(p::OptParams, lc::LoadCases, state)
    _design_loads(p) || return lc.F
    q = _frame_fef_stacked(p.model.cache.fef, state)       # 12·n_loads
    return lc.P .- lc.static .- p.model.cache.fef.scatter * (q .* lc.M)
end

# every frame element load's global fixed-end 12-vector at the state, stacked
# in cache.fef order (the unscattered form of Asap.assemble_Pf's frame part)
function _frame_fef_stacked(t, state)
    Q = map(eachindex(t.loads)) do j
        e = t.elem[j]
        ends = state.ends === nothing ? t.ends[j] : state.ends[e]
        Asap.element_fef_global(t.loads[j], state.sections[e], ends,
            _col3(state.X, t.i1[j]), _col3(state.X, t.i2[j]), t.rollangle[j])
    end
    return reduce(vcat, Q)
end

"""
    compliance(res::OptResults, p::OptParams, lc::LoadCases) -> Vector

Per-case external work `Fₖᵀuₖ` of a [`solve_cases`](@ref) result.
"""
compliance(res::OptResults, p::OptParams, lc::LoadCases) = vec(sum(res.F .* res.U; dims=1))

"""
    axial_force(res::OptResults, p::OptParams, lc::LoadCases) -> Matrix

Per-element, per-case axial force (n_el × n_cases), tension-positive — the
multi-case form of [`axial_force`](@ref)`(res, p)`, including each case's own
element-load axial fixed-end force `q̃₇`.
"""
function axial_force(res::OptResults, p::OptParams, lc::LoadCases)
    ΔX = res.X * transpose(p.Cinc)
    g = vec(res.EA ./ (res.L .^ 2))
    N = reduce(hcat, map(1:length(lc)) do k
        ΔU = reshape(res.U[:, k], 6, :)[1:3, :] * transpose(p.Cinc)
        vec(sum(ΔX .* ΔU; dims=1)) .* g
    end)
    size(p.Eload, 2) == 0 && return N
    _design_loads(p) || return N .+ lc.qa0
    Q = Asap.fef_local(p.model.cache, ModelState{eltype(res.X)}(res.X, res.sections, res.EA, res.ends))
    return N .+ p.Eload * (Q[7, :] .* lc.E)
end

"""
    axial_stress(res::OptResults, p::OptParams, lc::LoadCases) -> Matrix

Per-element, per-case axial stress (n_el × n_cases).
"""
axial_stress(res::OptResults, p::OptParams, lc::LoadCases) = axial_force(res, p, lc) ./ res.A

function Base.show(io::IO, ::MIME"text/plain", lc::LoadCases)
    print(io, "LoadCases  ($(length(lc)) case$(length(lc) == 1 ? "" : "s"): ",
          join((":" * String(c) for c in lc.names), ", "), ")")
end
