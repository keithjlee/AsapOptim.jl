# S4.2 gb_mma: publication stack vs Asap/AsapOptim v1.x

**Date**: 2026-07-17 · Apple Silicon (Keith's MacBook), Julia 1.11.9 for BOTH
stacks · scripts: `docs/benchmarks/gb_mma_old.jl` / `gb_mma_new.jl`

## Question

The DiffAnalysis_2024 paper's gradient-based benchmark
(`S4.2_minvolume_spaceframe/gb_mma.jl`) minimizes spaceframe volume under
displacement + stress constraints with MMA, stopping at a 300 s time limit.
How does the modern stack — Asap v1.1 + AsapOptim v1.0, NLopt driven
directly, reverse-mode objective + prepared ForwardDiff constraint
Jacobian + `CachedSolver` — compare on the identical problem?

(The paper environment was originally Julia 1.10.2; its Manifest has since
been re-locked on 1.11.3, and both stacks here run on 1.11.9 — apples to
apples on this machine, which is faster than the paper's: the old stack no
longer needs the full 300 s at the paper's `ftol_rel = 1e-3`.)

## Problem (identical on both sides, verified)

- 145 nodes / 512 elements / 339 free DOFs (the pinned publication fixture;
  displacement parity vs the publication solve: 1.1e-13)
- **540 design variables, mixed**: 28 independent top-node Z positions (with
  diagonal-mirror `CoupledVariable` pairs) + 512 element areas
- **596 constraint rows**: 145 one-sided vertical displacements + 451
  stresses (elements not spanning two supports), exactly the paper's
  unnormalized rows
- MMA (`LD_MMA`), `maxeval = 1000`, `maxtime = 300`
- Parity checks: `x_init` and bounds **bit-identical** between stacks;
  initial volume identical to the last digit (160.14539078784395); initial
  design feasible on both.

Old stack = the publication's own environment and driver
(`constrained_optimization`, Nonconvex + Zygote objective AND Zygote
constraint Jacobian). New stack = NLopt C API directly, Zygote objective
gradient, prepared ForwardDiff constraint Jacobian (596×540), one
`Asap.CachedSolver` shared by all passes.

## Results

| Run | Stack | Stop | Final volume | Wall time | Iterations | ms/iter |
|---|---|---|---|---|---|---|
| paper settings (`ftol_rel=1e-3`) | publication | FTOL | 2.6954 | 192.4 s | 347 | 554 |
| | **v1.x** | FTOL | 2.7370 | **48.5 s** | 235 | 206 |
| full budget (`ftol_rel≈0`) | publication | MAXTIME (300 s) | 2.6585 | 302.5 s | 548 | 552 |
| | **v1.x** | **MAXEVAL (1000)** | **2.6562** | **103.6 s** | 1000 | **104** |

All four runs end feasible (max constraint ≤ −4e-8). Volumes are local
MMA outcomes along different trajectories; the meaningful comparisons:

- **Per-iteration cost: 5.3× faster** (104 vs 552 ms at matched full-length
  runs). The constraint Jacobian dominates an MMA iteration, and forward
  mode + the shared factorization is where the time goes away (44,196
  factorization cache hits vs 1,851 refactorizations over the run).
- **Under the paper's stopping criteria, the new stack is no longer
  time-limited**: it exhausts the paper's own 1000-iteration cap in 103.6 s
  — a third of the budget — and its result (2.6562) is slightly better
  than what the publication stack reaches using the entire 300 s (2.6585).
  Raising `maxeval` would let it keep converging inside the same budget.
- At the paper's loose `ftol_rel = 1e-3` both stacks stop early at
  statistically similar quality (2.70 vs 2.74, 1.5% apart — different
  trajectories through a nonconvex landscape); the new stack gets there
  4.0× sooner.

## Why the trajectories differ

Formulation, settings, start point, and constraint rows are identical; the
optimizer differs only in how derivatives are computed (bit-differences at
the 1e-13 level compound over hundreds of nonconvex MMA iterations) and in
the driver layer (Nonconvex's wrapper vs the NLopt C API directly). Both
end feasible at comparable optima — the trajectories are siblings, not
copies.

## Reproduce

```bash
OUTDIR=/tmp/gbmma julia +1.11 docs/benchmarks/gb_mma_old.jl   # also dumps problem structure
OUTDIR=/tmp/gbmma julia +1.11 docs/benchmarks/gb_mma_new.jl
# full-budget variants:
FTOL_REL=1e-12 OUTDIR=... julia +1.11 docs/benchmarks/gb_mma_{old,new}.jl
```

The old script requires the publication repository at its OneDrive path
(see the script header); per-iteration objective/time histories are written
to `OUTDIR`.
