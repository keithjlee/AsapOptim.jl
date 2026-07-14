# AsapOptim.jl

Structural optimization layer over Asap.jl (`../Asap`). **v1.0 (branch `asap-v1`): rewritten on Asap's differentiable core.**

The package is now a thin, Zygote-FREE design-vector layer (~600 lines): variables (`SpatialVariable`/`AreaVariable`/`CoupledVariable`) compile into an `OptParams` (sparse scatter maps + the processed model); `solve_structure(x, p)` (aliases `solve_truss`/`solve_frame` — the core is unified) evaluates a design through `Asap.ModelState` → `Asap.solve` (pure, AD-transparent); `axial_force`/`axial_stress`/`compliance` extend Asap's accessors on `OptResults`; `GeometricProperties(x, p)` is the no-solve geometry path; `updatemodel(p, x)` materializes a design back into a solved `Model`.

There are NO custom rrules here — Asap's `AsapChainRulesExt` (activated by loading Zygote/ChainRulesCore) covers everything. The legacy v0.1.x implementation (functional re-assembly, all_inz sparsity hack, hand-written adjoints) is preserved unloaded in `legacy_v0/`.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'   # incl. Zygote-vs-FiniteDifferences gradient checks
```

## Known gaps (deliberate, documented)

- `SectionVariable` (full section parameterization) not yet ported; `RigiditySection`-based parameterization is the intended v1.0 idiom.

The force-density Network path IS ported: `QVariable` → `NetworkOptParams` → `solve_network` (forward parity with Asap's FDM solver verified at 1e-10; gradients via Asap's multi-RHS solve_free adjoint).

## Downstream consumers

`../DemandTransport2` (active research) consumes `solve_truss`, `TrussOptParams`, variables, `res.U`/`res.L`, `axial_force(res, p)`, `GeometricProperties`, `updatemodel`. Its port (Phase 5d): `TrussOptParams` still works (alias); `res.U` is now the FULL 6-DOF/node vector (vertical displacements at `U[3:6:end]`·pattern, not `[2:3:end]`).
