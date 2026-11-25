"""
    LoadCases

Several load cases of one optimization problem, solved together against ONE
stiffness factorization per design (a multi-right-hand-side solve). Built
from an [`OptParams`](@ref) and the case tags of the model's loads:

    lc = LoadCases(p)                       # every case in the model, first-appearance order
    lc = LoadCases(p; cases = [:dead, :snow])

A thin wrapper over Asap's `CaseLoads` (`Asap.case_loads`), which holds every
case's loads as plain data. Everything here is constant per optimization; the
design-dependent part of each case's load vector (`p.loads = :design`: frame
line and point loads and self-weight, quadrilateral and triangle shell
surface loads) is re-evaluated per design by [`solve_cases`](@ref).

A load combination is a linear combination of cases: solve the cases once
and combine `res.U * λ` (or build the combination as its own case by tagging
factored loads with one case name). `Asap.combination_weights(lc.cl, combo)`
turns an `Asap.LoadCombination` into λ.

# Fields
- `names::Vector{Symbol}`: the case tags, one column each below
- `cl::Asap.CaseLoads{Float64}`: the case table (nodal loads, frozen load
  vectors, every case's element loads)
- `P::Matrix{Float64}`: full-space nodal loads per case (n_dof × n_cases)
- `F::Matrix{Float64}`: full-space `P − Pf` per case at the BASE design (the
  frozen load vectors)
- `static::Matrix{Float64}`: design-independent fixed-end forces per case
  (element loads on `VariableElement`s)
- `E::Matrix{Float64}`: n_frame_loads × n_cases case membership of each frame
  element load, in `cl.fef` order
- `Eel::SparseMatrixCSC{Float64,Int}`: n_el × n_frame_loads, 1 where frame
  load `j` (`cl.fef` order) sits on element `e`
- `qa0::Matrix{Float64}`: per-element axial fixed-end force `q̃₇` per case at
  the base design (frozen loads)
"""
struct LoadCases
    names::Vector{Symbol}
    cl::Asap.CaseLoads{Float64}
    P::Matrix{Float64}
    F::Matrix{Float64}
    static::Matrix{Float64}
    E::Matrix{Float64}
    Eel::SparseMatrixCSC{Float64,Int}
    qa0::Matrix{Float64}
end

function LoadCases(p::OptParams; cases::AbstractVector{Symbol}=Asap.load_cases(p.model))
    isempty(cases) && throw(ArgumentError("no load cases"))
    model = p.model
    cl = Asap.case_loads(model; cases)          # validates the names; leaves the cache as it was
    nl = length(cl.fef.loads)
    nc = length(cl.names)
    E = zeros(nl, nc)
    for (j, k) in enumerate(cl.fcase)
        E[j, k] = 1.0
    end
    Eel = sparse(cl.fef.elem, collect(1:nl), ones(nl), length(model.elements), nl)
    # frozen axial fixed-end forces per element and case, at the model's own state
    Q0 = Asap.fef_local(cl, Asap.extract_state(model))
    qa0 = nl == 0 ? zeros(length(model.elements), nc) : Matrix(Eel * (Q0[7, :] .* E))
    return LoadCases(cl.names, cl, cl.P, cl.F0, cl.static, E, Eel, qa0)
end

Base.length(lc::LoadCases) = length(lc.names)

"""
    solve_cases(x, p::OptParams, lc::LoadCases) -> OptResults

Evaluate a design for every case in `lc` with one stiffness assembly and one
factorization: like [`solve_structure`](@ref), but `res.U` and `res.F` are
matrices (n_dof × n_cases, column `k` = case `lc.names[k]`). Positions,
sections and lengths are shared. Differentiable through the same paths
(ForwardDiff Duals, ChainRules engines — Asap's multi-RHS `solve_free`
rules and its case-load rules).

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
    F = Asap._case_F(lc.cl, state, p.loads)        # frozen, or P − Pf(x) per case
    cm = cache.constraint
    U = cm === nothing ? cache.free_embed * Asap._pure_solve(K, F[free, :], p.solver) :
        cache.free_embed * Asap._constrained_solve(cm, K, F[free, :], p.solver)
    return OptResults(U, X, A, _element_lengths(X, p), sections, EAvec, F, ends)
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
    size(lc.E, 1) == 0 && return N
    p.loads === :design || return N .+ lc.qa0
    Q = Asap.fef_local(lc.cl, ModelState{eltype(res.X)}(res.X, res.sections, res.EA, res.ends))
    return N .+ lc.Eel * (Q[7, :] .* lc.E)
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
