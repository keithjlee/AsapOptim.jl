# AD Verification & Performance: Asap v1.0 vs the Publication Stack

**Date**: 2026-07-13 · **Machine**: Apple Silicon (Keith's MacBook), Julia 1.11 for BOTH stacks

## Question

Does automatic differentiation work — provably correctly — through the new
Asap v1.0 / AsapOptim v1.0 stack, and how do speed and memory compare to the
legacy stack (registered Asap 0.2.1 + the publication's differentiable layer,
i.e. the DiffAnalysis_2024 vendored AsapOptim) on a real problem?

## Setup

The **S4.2 minimum-volume spaceframe** from the DiffAnalysis_2024 paper
repository: **145 nodes, 512 truss elements, 339 free DOFs**.

- The legacy stack builds the model **verbatim through the publication's own
  `init_problem.jl`** in its own pinned environment.
- The new stack rebuilds the *identical* structure from the serialized
  publication fixture (`Asap/test/characterization/fixtures_diffanalysis.jl`) —
  verified: displacements match the publication solve to **1.1 × 10⁻¹³**.
- Same design parameterization on both sides: one absolute area variable per
  element (512 design variables), identical objectives:
  - compliance `Fᵀu` (through the solve)
  - volume `Σ AᵢLᵢ` (geometry only)
- Benchmarks are medians of 10 (solves) / 5 (gradients) runs after warm-up,
  with GC between runs. Scripts in `docs/benchmarks/` (`bench_old.jl`,
  `bench_new.jl`) reproduce everything.

## AD correctness — verified at four levels

| Check | Result |
|---|---|
| Zygote vs central finite differences (5-point), directional derivative of **compliance** | rel. err **1.9 × 10⁻⁹** |
| Zygote vs finite differences, **volume** | rel. err **8.6 × 10⁻¹⁵** |
| **New-stack gradient vs legacy-stack gradient**, compliance (512 components) | max rel. diff **3.6 × 10⁻¹³**, cosine 1.0 |
| New vs legacy gradient, volume | **bit-identical** |

Plus the standing test suites: Asap core (gradients w.r.t. node positions,
areas, and semi-rigid connection stiffness vs finite differences), AsapOptim
(coupled variables, grouped areas, mixed frame+truss, FDM network path).

**Conclusion: yes — AD through the new stack is correct**, and reproduces the
publication-era gradients on the paper's own structure to machine precision.

## Performance

| Metric (median) | Legacy (0.2.1 + paper layer) | **New v1.0** | Ratio |
|---|---|---|---|
| `solve!` (assembly + factorize + solve + post-process) | 1.03 ms / 4.30 MB | **0.88 ms** / 4.36 MB | 1.2× faster |
| stiffness assembly alone | 0.199 ms / 1.56 MB | **0.130 ms** / 1.37 MB | 1.5× faster |
| differentiable forward solve | 0.837 ms / 2.32 MB | **0.413 ms** / 1.58 MB | 2.0× faster |
| ∇ compliance (Zygote, 512 vars) | 1.41 ms / 7.60 MB | **2.77 ms** / 12.2 MB | 2.0× slower |
| ∇ volume | 0.27 ms / 0.11 MB | **2.90 ms** / 18.3 MB | ~10× slower |

Context for the two gradient rows: the legacy layer achieved its speed with
**hand-derived analytic rrules for every operation** (fixed derivative
matrices, hard-coded to the 3-DOF truss layout). The new stack differentiates
the *generic* pipeline — the same code path that handles frames, semi-rigid
joints, springs, and super-elements — with exactly one performance rule
(analytic `truss_stiffness` pullback). Being within 2× of a fully
hand-differentiated implementation while keeping full generality is the
intended trade; both gradients run at hundreds of evaluations/second at this
problem size.

## The optimization journey (engineering notes for future work)

The first working gradient took **6,631 ms and churned 612 GB**. Three
structural fixes brought it to 2.8 ms (2,400×). These patterns matter for any
future differentiable code in this ecosystem:

1. **Never read mutable-struct fields inside a differentiated closure.**
   `el.nodeStart.index` per element made Zygote build and accumulate tangent
   structures for the entire model object graph per access. Fix: plain-data
   mirrors (`i1`, `i2`, `Ψs`, `endss` on `ElementGroup`; `i1`, `i2`,
   `base_sections`, `Evec` on `OptParams`). → 6,631 ms → 102 ms.
2. **Avoid per-element closures and generic `getindex` pullbacks; batch.**
   Truss groups now assemble with one incidence matmul (`X * Cincᵀ`) and 36
   broadcasted row expressions — a handful of dense array ops instead of
   thousands of scalar graph nodes. → 102 → 50 ms (with the analytic
   `truss_stiffness` rule for the remaining per-element paths).
3. **Struct getfield chains per element are as bad as struct reads.**
   `map(i -> EA(sections[i]), indices)` alone cost 42 ms. Fix: carry `EA` as
   a plain precomputed vector on `ModelState`/`OptResults`. → 50 → **2.8 ms**.

### Update 2026-07-17: PARITY REACHED AND EXCEEDED — this supersedes the headline table

After the masked-scatter fix (below), the incidence-matmul geometry
batching, and Asap v1.1.1, the original "2× / 10× slower" verdict is
obsolete. Same scripts, same machine, Julia 1.11, fresh same-session runs
of BOTH stacks (legacy = the publication environment, as in the original
table):

| Metric (median) | Legacy (0.2.1 + paper layer) | **New v1.x** | Ratio |
|---|---|---|---|
| `solve!` | 0.95 ms | **0.58 ms** | 1.6× faster |
| stiffness assembly alone | 0.186 ms | **0.045 ms** | 4.1× faster |
| differentiable forward solve | 0.78 ms | **0.40 ms** | 1.9× faster |
| ∇ compliance (Zygote, 512 vars) | 1.29 ms | 1.53 ms | 1.19× slower |
| ∇ compliance (**Enzyme**, 512 vars) | — | **0.83 ms** | **1.55× faster** |
| ∇ volume | 0.229 ms | **0.158 ms** | 1.45× faster |

And on the workloads that dominate real constrained optimization
(measured in `docs/AD_BACKENDS.md` / `docs/BENCHMARK_GB_MMA.md`):

- **512×512 stress-constraint Jacobian**: legacy ~400 ms → implicit
  differentiation **3.8 ms** (~105× faster).
- **The publication's own gb_mma benchmark end-to-end**: 552 → 57
  ms/iteration; its 1000-iteration budget completes in 56.8 s vs the
  publication stack spending 300 s on 548 iterations.

The sole remaining legacy win is Zygote-driven ∇compliance at +19% — the
documented price of the generic pipeline under the interpreter-based
engine, erased by switching to Enzyme (and irrelevant to constrained
optimization, where Jacobians dominate). AD checks unchanged
(1.9e-9 / 8.6e-15).

Status of the "known further opportunities" listed at the bottom:
lengths batching DONE (the ∇volume flip); the `solve_free` cotangent
projection DONE; Enzyme/Mooncake DONE (see AD_BACKENDS.md); batched
frame-element assembly remains open.

### Update 2026-07-16: masked-scatter broadcast fix

The per-element mask comprehensions in `_design_state`,
`GeometricProperties`, and `solve_network`
(`[mask[i] ? var[i] : base[i] for i ...]`) were pattern-1/2 offenders in
disguise: each indexed read costs a one-hot `getindex` pullback (~35×
a broadcasted `ifelse` in isolation). Replaced with
`ifelse.(mask, Sa * x, A0)`; sections also now skip the per-element
rebuild when no area variables exist. Same machine, Julia 1.11, same
fixture — gradients bit-identical, AD checks unchanged:

| Metric (median) | before | after | vs legacy 1.41 / 0.27 ms |
|---|---|---|---|
| ∇ compliance (512 area vars) | 2.77 ms | **1.78 ms** | 1.26× |
| ∇ volume | 2.90 ms | **1.61 ms** | 6× |

On the small 24-variable example-1 truss the effect is larger (fixed
overheads dominate there): `Zygote.withgradient` of compliance went
945 → **373 μs** (legacy stack: 161 μs), ∇ volume 627 → **72 μs**.

### Known further opportunities (not yet done)

- `∇ volume` (now 1.61 ms) is dominated by the per-element
  `_element_lengths` closure — batchable with the same incidence-matmul
  trick (expect ~0.3 ms).
- Frame-element groups still use the per-element kernel path; a batched
  frame assembly (or an analytic `frame_stiffness` pullback) would matter for
  frame-dominated optimization.
- The `sparse_from_pattern`/`solve_free` pullback materializes the dense
  `−λuᵀ` outer product; a lazy/thunked projection straight onto the sparsity
  pattern would cut gradient allocations further.
- Enzyme/Mooncake on the pure path (the mutation-free design was built for
  this) — untested.

## Reproduce

```bash
julia +1.11 docs/benchmarks/bench_new.jl    # new stack (temp env, dev'd local packages)
julia +1.11 docs/benchmarks/bench_old.jl    # publication stack (its own pinned env)
```
