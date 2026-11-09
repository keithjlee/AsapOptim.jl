# Examples

Runnable scripts in the repository's [`examples/`](https://github.com/keithjlee/AsapOptim.jl/tree/main/examples) folder:

| File | What it shows |
|---|---|
| `examples/truss-optimization1.jl` | geometry optimization with a compliance objective (the [Quick start](quickstart.md)) |
| `examples/truss-optimization2.jl` | constrained minimum volume: geometry and sizing, with MMA |
| `examples/truss-optimization3-nlopt.jl` | the recommended pattern: NLopt driven directly, a reverse-mode objective, a forward-mode constraint Jacobian and `CachedSolver` (about 9× faster per run than the all-Zygote formulation) |
| `examples/vierendeel-geometry.jl` | frame optimization: a Vierendeel truss profile, volume-scaled compliance, hand-rolled projected gradient descent |
| `examples/joint-stiffness.jl` | [`JointVariable`](@ref): the cheapest distribution of semi-rigid connection stiffness that meets a drift limit |
| `examples/ad_backends/` | the same problem differentiated with Zygote, Mooncake and Enzyme |

Engineering reports with full benchmark tables are in the repository's [`docs/`](https://github.com/keithjlee/AsapOptim.jl/tree/main/docs) folder: AD backend comparisons, AD verification against the published implementation, and the publication `gb_mma` benchmark.
