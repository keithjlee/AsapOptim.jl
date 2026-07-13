# AsapOptim.jl

Differentiable structural optimization layer over Asap.jl (`../Asap`). Zygote + ChainRulesCore + LinearSolve. Pinned to `Asap = "0.2"`.

**Asap is undergoing a v1.0 modernization (`../Asap/docs/MODERNIZATION.md`); this package migrates in lockstep at Phase 5a.** The migration *deletes* most of this package's Functions layer: Asap's new AD-first core provides the pure functional assembly/solve path natively, so `Functions/{K,Kframe,Ktruss,Rframe,Rtruss}.jl`, the `all_inz` sparsity hack in `Types/Utilities.jl`, and most hand-written rrules go away. What survives: `Types/Variables/`, `Types/Indexers/`, `Constraints.jl`, objective functions, and only rrules the core doesn't cover.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Architecture

- `Types/Variables/` — design variables (`SpatialVariable`, `AreaVariable`, `SectionVariable`, `QVariable`, `CoupledVariable`)
- `Types/Indexers/` — map variable positions to model-field indices
- `Types/Parameters/` — `TrussOptParams`/`FrameOptParams`/`NetworkOptParams`: differentiable snapshots of a processed Asap model (extract geometry, section values, DOF maps, and the sparsity pattern of `model.S`)
- `Functions/` — pure re-implementations of Asap's pipeline (geometry → transformations → element k → global K → solve) with hand-written `ChainRulesCore.rrule`s; `_noadjoint` twins exist for finite-difference validation
- `Functions/Solve.jl` — `solve_truss`/`solve_frame` entry points

## Critical coupling to Asap (breaks if Asap changes)

- `all_inz` / `get_inz` (`Types/Utilities.jl`) reverse-engineer `model.S.colptr/rowval/nzval` with hardcoded 3 (truss) / 6 (frame) DOFs per node.
- `Kframe.jl:127-198` / `Ktruss.jl` are **verbatim numeric copies** of Asap's `local_K` family — must stay bit-identical. (`FreeFree` release has no `k` method here — latent bug.)
- Hardcoded force/DOF indices: truss axial = `fvecs[2]`, frame axial = `fvecs[7]`, truss stress rows `4:6`.
- `updatemodel`/`updatenetwork` (`PostProcessing.jl`) rebuild Asap objects via positional constructors.

## Downstream consumers

`../DemandTransport2` (active research) consumes `solve_truss`/`solve_network`, `TrussOptParams`/`NetworkOptParams`, `SpatialVariable`/`QVariable`/`CoupledVariable`, results fields `res.U/.L/.X/.Y/.Z/.Q`, `axial_force(res, p)`, `GeometricProperties`, and `updatemodel` — inside Zygote-differentiated objectives. The Phase 5a rework must keep these (or equivalents) working and differentiable w.r.t. node positions and force densities.

## Known bugs

- `SectionVariable.iglobal::Float64` (should be Int).
- `f_axial_new` pullback has an early `return` inside its per-element loop — only one element's gradient accumulates.
