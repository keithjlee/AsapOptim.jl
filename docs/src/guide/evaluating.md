# Evaluating designs

## `solve_structure(x, params) -> OptResults`

[`solve_structure`](@ref) is the full differentiable analysis. It scatters the design vector into positions, sections and joint stiffnesses, builds the load vector, runs Asap's pure solve, and packages the results. `solve_truss` and `solve_frame` are aliases: the core is unified, so one function serves trusses, frames and mixed models.

[`OptResults`](@ref) fields are all plain data, so any scalar function of them is differentiable back to `x`:

- `U`: the global displacement vector, **6 DOFs per node** (`ux, uy, uz, θx, θy, θz`). Node `i`'s translations are `U[6(i-1) .+ (1:3)]`, all y-displacements are `U[2:6:end]`, and all z-displacements are `U[3:6:end]`.
- `X`: evaluated node positions (3 × n_nodes).
- `A`, `L`: per-element areas and lengths.
- `sections`, `EA`: evaluated sections and axial rigidities.
- `F`: the load vector the design was solved for. It varies with the design when members carry element loads.
- `ends`: evaluated end conditions (joint variables), or `nothing`.

Built-in reductions of an `OptResults`:

- `compliance(res, params)`: external work `Fᵀu`, the canonical smooth stiffness objective.
- `axial_force(res, params)`: per-element axial force, tension-positive, including the axial component of element loads.
- [`axial_stress`](@ref)`(res, params)`: axial force over the evaluated area.

## Several load cases: `LoadCases` and `solve_cases` (AsapOptim ≥ 1.4)

Design problems usually check several load cases (dead, snow, wind, …) at the same design. Instead of one `OptParams` per case, collect the cases once and solve them together against **one** stiffness assembly and factorization (a multi-right-hand-side solve):

```julia
params = OptParams(model, variables)          # loads tagged with their case: NodeForce(…; case = :snow)
lc  = LoadCases(params)                        # every case in the model, in first-appearance order
lc  = LoadCases(params; cases = [:dead, :snow, :wind])

res = solve_cases(x, params, lc)               # res.U, res.F: n_dof × n_cases (column k = lc.names[k])
c   = compliance(res, params, lc)              # per-case Fₖᵀuₖ
N   = axial_force(res, params, lc)             # n_el × n_cases, tension-positive
σ   = axial_stress(res, params, lc)

# a load combination is a linear combination of the cases
λ   = [1.0, 1.0, 0.0]                          # D + S
Nc  = N * λ
Uc  = res.U * λ
```

[`LoadCases`](@ref) wraps Asap's `CaseLoads` (Asap ≥ 1.5): everything that is constant during the optimization, i.e. the per-case nodal loads, the frozen load vectors, static fixed-end forces and every case's element loads. [`solve_cases`](@ref) re-evaluates the design-dependent part (`loads = :design`: frame line and point loads and self-weight, quadrilateral and triangle shell surface loads) per case at each design, so `axial_force` includes each case's own element-load axial component. A `LoadCombination` becomes column weights with `Asap.combination_weights(lc.cl, combo)`.

`solve_cases` is differentiable through the same paths as `solve_structure` (ForwardDiff, Zygote and other ChainRules engines). Its results equal one `OptParams` per case, and its gradients equal the sum of the single-case gradients; the test suite checks both.

## `GeometricProperties(x, params)`

[`GeometricProperties`](@ref) is a geometry-only evaluation, returning positions `X`, lengths `L` and areas `A` with **no structural solve**. Use it for objectives like material volume (`dot(geo.A, geo.L)`), where paying for a solve and its gradient would be waste. For shells, [`shell_volume`](@ref)`(x, params)` gives Σ t·area at the design, also without a solve.

## `updatemodel(params, x) -> Model`

[`updatemodel`](@ref) materializes a design back into the reference `Asap.Model`. It writes positions, sections and joint stiffnesses, re-solves, and returns the ordinary mutable model, ready for post-processing, force recovery, visualization or export.
