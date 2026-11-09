# Migrating from v0.1

- `TrussOptParams` and `FrameOptParams` still work. They are aliases of the unified [`OptParams`](@ref), and `solve_truss`/`solve_frame` alias [`solve_structure`](@ref).
- `res.U` is the **full 6-DOF-per-node** vector for every model (v0.1's truss path was 3-DOF): vertical displacements are `U[2:6:end]`, not `U[2:3:end]`.
- Results carry positions as a 3 × n matrix `res.X` instead of separate `X`/`Y`/`Z` vectors.
- The load vector of an evaluated design is `res.F` (with element loads it varies with the design); `params.F` is the base-design vector. Or use `compliance(res, params)`.
- Custom rrules, the functional re-assembly path and the truss/frame code split are gone: differentiability comes from Asap's core. The old implementation is preserved, unloaded, in `legacy_v0/` for reference.
- Section parameterization is available through [`SectionVariable`](@ref) (`:A`, `:Ix`, `:Iy`, `:J`).
