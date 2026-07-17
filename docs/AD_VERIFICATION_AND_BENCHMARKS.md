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
