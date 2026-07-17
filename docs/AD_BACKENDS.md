# AD Backend Support: Zygote, Mooncake, Enzyme

**Date**: 2026-07-13 · Julia 1.11, Apple Silicon · via DifferentiationInterface.jl

## Question

Does the new Asap + AsapOptim work with **Mooncake** and **Enzyme** — the
actively developed AD engines — so the suite can escape its Zygote
dependency? (The legacy AsapOptim worked with neither.)

## Answer

**Yes — all three backends now work, out of the box, to identical
accuracy.** Two new package extensions were required (loaded automatically;
no user code changes beyond backend selection):

- `AsapMooncakeExt` — bridges Asap's ChainRules rules to Mooncake
  (`Mooncake.@from_rrule`), supplies the sparse-cotangent adapter Mooncake's
  interop lacks, and declares CHOLMOD factorizations tangent-free.
- `AsapEnzymeExt` — imports the two rules Enzyme cannot traverse (the
  frozen-pattern sparse construction and the CHOLMOD solve) via
  `Enzyme.@import_rrule`, registered at load time (`__init__`) because
  Enzyme's own ChainRules extension isn't active during precompile.

Asap's core rules were also made **representation-flexible**: the
`solve_free` K-cotangent is now returned projected onto the frozen sparsity
pattern (no dense n×n outer product — a win for every backend), and
`sparse_from_pattern`'s pullback accepts both matrix and structural-tangent
cotangents.

## Usage

```julia
using DifferentiationInterface
using Asap, AsapOptim

obj(x) = compliance(solve_structure(x, p), p)

# Zygote — works as before
import Zygote
gradient(obj, AutoZygote(), x0)

# Mooncake — nothing special needed
import Mooncake
gradient(obj, AutoMooncake(config = nothing), x0)

# Enzyme — two annotations: runtime activity (the model is captured constant
# state) and Const function treatment (closure captures are data)
import Enzyme
backend = AutoEnzyme(
    mode = Enzyme.set_runtime_activity(Enzyme.Reverse),
    function_annotation = Enzyme.Const)
gradient(obj, backend, x0)
```

## Correctness

Compliance gradient, max relative deviation from 5-point central finite
differences:

| Problem | Zygote | Mooncake | Enzyme |
|---|---|---|---|
| small truss (3 mixed vars) | 8.7e-13 | 8.7e-13 | 8.7e-13 |
| publication S4.2 spaceframe (512 area vars) | 1.30e-9 | 1.30e-9 | 1.30e-9 |

All three agree with each other to machine precision.

## Performance (median gradient evaluation, warm)

| Problem | Zygote | Mooncake | **Enzyme** | legacy stack* |
|---|---|---|---|---|
| small (3 vars) | 0.76 ms / 0.40 MB | 0.32 ms / 0.08 MB | **0.19 ms / 0.11 MB** | — |
| spaceframe (512 vars) | 3.15 ms / 11.8 MB | 4.47 ms / 3.9 MB | **1.04 ms / 3.6 MB** | 1.41 ms / 7.6 MB |

\* legacy = Asap 0.2.1 + the publication's fully hand-differentiated layer
(see `AD_VERIFICATION_AND_BENCHMARKS.md`).

Headlines:

- **Enzyme is the fastest backend on every problem — and beats the legacy
  hand-written-adjoint stack** (1.04 vs 1.41 ms) while running the fully
  generic pipeline. The escape from Zygote comes with a speedup, not a tax.
- Mooncake matches Zygote's accuracy with ~3× less memory at scale; it wins
  on small problems and currently trails Zygote in time on larger ones.
- Zygote continues to work unchanged.

## Update 2026-07-16: example-1 truss, backend vs legacy stack

After the masked-scatter broadcast fix (AsapOptim) and the `truss_stiffness`
static-block fix (Asap), the line-95 compliance gradient of
`examples/truss-optimization1.jl` (24 spatial variables, 47 elements),
`@benchmark` medians, all gradients agreeing to 1.7e-10 vs finite
differences:

| Stack / backend | median | memory |
|---|---|---|
| legacy (Asap 0.2.2 + AsapOptim 0.1.3), Zygote | 179 μs | 0.49 MB |
| v1.0, Zygote | 446 μs | 0.78 MB |
| v1.0, Mooncake (DI, prepared) | 253 μs | 0.29 MB |
| **v1.0, Enzyme (DI, prepared)** | **189 μs** | 0.39 MB |

Enzyme on the generic v1.0 path matches the fully hand-differentiated
legacy stack (within 6%) even on this small, truss-only,
legacy-favorable problem — consistent with the spaceframe result above,
where Enzyme beats legacy outright.

And the S4.2 spaceframe (512 area variables), same conditions — the
legacy baseline here is the REGISTERED old stack (Asap 0.2.2 +
AsapOptim 0.1.3 with current Zygote), which is faster than the
publication's pinned environment measured in the original table above
(0.83 vs 1.41 ms), so this is the tougher comparison:

| Stack / backend | Julia 1.12.6 | Julia 1.11 |
|---|---|---|
| legacy 0.1.3, Zygote | 0.832 ms / 7.5 MB | — |
| v1.0, Zygote | 1.433 ms / 8.5 MB | 1.43 ms / 8.7 MB |
| v1.0, Mooncake (DI, prepared) | 4.35 ms / 3.9 MB | 4.46 ms / 3.9 MB |
| v1.0, Enzyme (DI, prepared) | 1.69 ms / 4.3 MB | **0.832 ms** / 3.5 MB |

All backends agree with the finite-difference directional derivative to
1.9e-9. Two observations:

- On Julia 1.11, Enzyme matches the legacy hand-differentiated stack
  EXACTLY (0.832 ms) at 512 variables — zero generality tax.
- **Enzyme is 2× slower on Julia 1.12 than 1.11 for identical code**
  (1.69 vs 0.832 ms; Zygote and Mooncake are version-stable). This is an
  Enzyme×Julia-1.12 codegen interaction, not an Asap/AsapOptim issue —
  worth re-checking as Enzyme releases catch up to 1.12. On 1.12,
  Zygote (1.43 ms) is currently the fastest backend at this scale.

## Update 2026-07-16: forward mode + the full backend × version × task matrix

With Asap v1.1 (solver seam), the geometry batching (`Cinc` incidence
matmuls for lengths AND `axial_force`), and the new forward-mode support
(`AsapForwardDiffExt` + `solve_free` frule), the full comparison —
`@benchmark` medians via DifferentiationInterface with prepared
gradients/Jacobians; legacy = registered Asap 0.2.2 + AsapOptim 0.1.3;
all AD results agree to ≤1e-13. Scripts: `matrix_new.jl`,
`matrix_legacy.jl`, `matrix_enzyme_forward_probe.jl`.

**Compliance GRADIENT (reverse mode's home turf):**

| Problem | backend | Julia 1.11 | Julia 1.12 |
|---|---|---|---|
| ex1 truss (24 spatial vars) | legacy-Zygote | 0.19 ms | 0.22 ms |
| | v1.0 Zygote | 0.46 | 0.45 |
| | v1.0 Mooncake | 0.25 | 0.26 |
| | v1.0 Enzyme-rev | **0.11** | 0.19 |
| | v1.0 ForwardDiff | 0.19 | 0.19 |
| spaceframe (512 area vars) | legacy-Zygote | 0.94 | 0.90 |
| | v1.0 Zygote | 1.35 | 1.34 |
| | v1.0 Mooncake | 5.02 | 4.41 |
| | v1.0 Enzyme-rev | **0.83** | 1.66 |
| | v1.0 ForwardDiff | 56.2 | 54.8 |

**Full stress-constraint JACOBIAN (n_el rows × n_vars cols):**

| Problem | backend | Julia 1.11 | Julia 1.12 | memory |
|---|---|---|---|---|
| ex1 (47 × 24) | legacy-Zygote | 8.4 ms | 8.6 ms | 19 MB |
| | v1.0 Zygote | 28.8 | 27.3 | 311 MB |
| | v1.0 Mooncake | 12.6 | 13.3 | 15 MB |
| | v1.0 ForwardDiff | **0.21** | **0.20** | 2 MB |
| spaceframe (512 × 512) | legacy-Zygote | 400 | 528 | 3.4 GB |
| | v1.0 Zygote | 7587 | 8108 | 614 GB churn |
| | v1.0 Mooncake | 2575 | 2487 | 2.1 GB |
| | v1.0 ForwardDiff | **58** | **61** | 0.5 GB |

Headlines:

- **Forward mode owns constraint Jacobians, decisively.** A forward
  Jacobian costs the same as a forward gradient — extra outputs are
  free — so ForwardDiff beats every reverse backend AND the legacy
  hand-differentiated stack: 40× vs legacy on the small problem, **7×
  vs legacy and 130× vs v1.0-Zygote at n_in = n_out = 512**. The
  intuition "structural optimization has more constraint outputs than
  design inputs, so forward should win" holds even at n_in = n_out.
- **Reverse mode still owns scalar objectives at scale** (spaceframe
  gradient: ForwardDiff 55 ms vs Zygote 1.3 ms). Rule of thumb:
  reverse (Enzyme-rev on 1.11 / Zygote on 1.12) for objectives,
  ForwardDiff for constraint Jacobians.
- Enzyme-rev on 1.11 now beats the legacy stack on the 512-var gradient
  (0.83 vs 0.94 ms); its 1.12 penalty (2×) persists.
- Enzyme-rev Jacobians FAIL via DifferentiationInterface ("Conversion
  of boxed type Vector{Float64}") — a DI↔Enzyme interop issue, not an
  Asap rule problem (gradients work).
- **Enzyme FORWARD mode is blocked upstream** on both Julia versions
  (0.13.186: compiler assertion `AdjointGenerator.h:318` on 1.12, a
  thrown compile error on 1.11). The `solve_free` frule import is in
  place and kernel-verified; re-probe on Enzyme releases with
  `matrix_enzyme_forward_probe.jl`.
- The geometry batching also moved reverse mode: spaceframe ∇volume
  1.61 → **0.19 ms**, ∇compliance 1.78 → 1.34 ms (Zygote, 1.12).

## Update 2026-07-17: Enzyme forward RESOLVED (native rule); solver seam; CachedSolver

Four changes landed together (Asap + AsapOptim; all suites green — 2718/57):

1. **Enzyme forward mode now works** (Julia 1.11, default solver). Root
   cause of the earlier failures was never Enzyme's forward core: the
   `@import_frule` bridge returned silently wrong tangents under runtime
   activity (a constant argument's shadow ALIASES its primal, which the
   generic bridge reads as a real tangent, `ΔF = F`), and the one BLAS
   call in the pipeline (`dot` in `compliance`) is unsupported in
   forward+runtime-activity mode. Fixes: a NATIVE `EnzymeRules.forward`
   rule in `AsapEnzymeExt` that detects shadow aliasing explicitly, and
   `sum(F .* u)` instead of `dot`. Verified to machine precision against
   Zygote (gradient and full Jacobian, batch widths 1–16).
2. **Solver selection reaches the differentiable path**:
   `OptParams(model, vars; solver = ...)` — `nothing` (CHOLMOD), any
   LinearSolve algorithm, or `Asap.CachedSolver()`. Values and gradients
   are backend-invariant (tested: KLU, KrylovJL_CG, CachedSolver, across
   Zygote AND ForwardDiff).
3. **CachedSolver** shares one factorization across every pass at a
   design iterate (obj gradient + constraint Jacobian + all ForwardDiff
   chunks). Spaceframe: Jacobian 54 → 43 ms, Zygote gradient 1.34 →
   1.06 ms. Factorization is only ~20% of these workloads at 339 DOFs —
   the cache's share grows with model size (Dual-K assembly per chunk
   dominates the rest; the next lever is a hand implicit-diff Jacobian
   for area problems, where K is linear in the design).
4. **NLopt-direct example** (`examples/truss-optimization3-nlopt.jl`):
   Zygote objective + prepared-ForwardDiff constraint Jacobian +
   CachedSolver → **500 MMA iterations in 6.6 s** where the
   Nonconvex/Zygote-Jacobian formulation needs its full 60 s budget.

Updated Jacobian standings (stress constraints, Julia 1.11, medians):

| Problem | ForwardDiff (+CachedSolver) | Enzyme-fwd (default solver) | Zygote |
|---|---|---|---|
| ex1 (47×24) | **0.16 ms** | 0.30 ms | 28.8 ms |
| spaceframe (512×512) | **43 ms** | 152 ms | 7,587 ms |

Guidance unchanged in direction, sharpened in detail: **ForwardDiff for
constraint Jacobians** (now also the fastest, not just the most robust),
reverse for scalar objectives, Enzyme-forward as an independent
correctness check. Enzyme limitations that remain upstream: 1.12 aborts
(compiler assertion) in forward mode, and Enzyme×CachedSolver aborts on
the mutable solver struct — with Enzyme use the default solver.

## Update 2026-07-17 (later): implicit-diff Jacobians — another 10×

`solution_tangents(x, p)` + `axial_force_jacobian`/`axial_stress_jacobian`
implement the implicit-function theorem directly: element-local
pseudo-loads (closed-form truss-area fast path; generic one-partial
ForwardDiff seed through the SAME kernels for everything else — spatial,
frame area, joints, and any future kernel-input variable), chained through
the compiled scatter maps (couplings/mixing structural), one factorization
+ one multi-RHS backsolve. Machine-precision parity (≤2e-14) with
full-pipeline ForwardDiff and Zygote on mixed problems including joints.

| Stress-constraint Jacobian (Julia 1.11) | time | memory |
|---|---|---|
| v1.0 Zygote | 7,587 ms | 614 GB churn |
| ForwardDiff (prepared, CachedSolver) | 39–43 ms | 0.5 GB |
| **implicit (`solution_tangents`)** | **3.9 ms** (512×512) / **4.5 ms** (gb_mma 596×540 mixed) | **26–34 MB** |

That is ~10× over prepared ForwardDiff and >100× over the LEGACY
hand-differentiated stack's 400 ms. Under NLopt-MMA the optimizer's
~150 ms/iter subproblem now dominates outright — constraint
aggregation/screening is the remaining lever. Plain-Float64 path for
optimizer callbacks (not itself AD-transparent); ForwardDiff/Zygote remain
the reference implementations in the test suite.

## Gotchas & their fixes (documented for future maintenance)

1. **Neither Mooncake nor Enzyme consumes ChainRules rules automatically** —
   both extensions exist precisely to register Asap's rules with each
   engine's native rule system.
2. **Mooncake + captured CHOLMOD factorizations**: an objective closing over
   a previously solved model captures a `CHOLMOD.Factor`; Mooncake tried to
   build tangent storage for it and failed internally. Fixed by declaring
   `Mooncake.tangent_type(::Type{<:CHOLMOD.Factor}) = NoTangent` (opaque
   foreign solver state is never differentiable).
3. **Enzyme + closures**: without `function_annotation = Enzyme.Const`,
   Enzyme treats the objective closure's captured `OptParams` as
   differentiable state and throws `EnzymeMutabilityException`. With the
   annotation (and runtime activity), everything works.
4. **`Enzyme.@import_rrule` in package extensions** must run in `__init__`
   (its implementation lives in Enzyme's own ChainRulesCore extension, which
   is not reliably loaded during a downstream extension's precompile).
5. A debugging red herring worth recording: a `ZeroPivotException` during
   backend testing turned out to be a **mechanism in the test model**
   (two collinear truss bars with a free transverse DOF), not an AD bug.
   Check `model.cache.partition` and equilibrium before blaming the engine.

## Not yet explored

- Enzyme forward mode / batched (vector) modes.
- Mooncake performance tuning (its `@from_rrule` bridge may add overhead vs
  native `rrule!!`s; a native rule for `stiffness_entries` would be the next
  step if Mooncake becomes the primary engine).
- Second derivatives / Hessians on any backend.

## Reproduce

`docs/benchmarks/adbench.jl` (temp environment; dev's local Asap/AsapOptim;
Julia 1.11).
