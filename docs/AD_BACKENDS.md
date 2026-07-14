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
