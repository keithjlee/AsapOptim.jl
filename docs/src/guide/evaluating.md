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

## `GeometricProperties(x, params)`

[`GeometricProperties`](@ref) is a geometry-only evaluation, returning positions `X`, lengths `L` and areas `A` with **no structural solve**. Use it for objectives like material volume (`dot(geo.A, geo.L)`), where paying for a solve and its gradient would be waste.

## `updatemodel(params, x) -> Model`

[`updatemodel`](@ref) materializes a design back into the reference `Asap.Model`. It writes positions, sections and joint stiffnesses, re-solves, and returns the ordinary mutable model, ready for post-processing, force recovery, visualization or export.
