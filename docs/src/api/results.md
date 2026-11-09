# Results and evaluation

Evaluating a design vector (`solve_structure`, `GeometricProperties`), the results it returns, the reductions over them (`compliance`, `axial_force`, `axial_stress`), and materializing a design back into a model (`updatemodel`).

```@autodocs
Modules = [AsapOptim]
Pages = ["results.jl", "update.jl"]
Private = false
```

AsapOptim extends two of Asap's accessors on `OptResults` (the functions belong to Asap, so they are listed here by method):

```@docs
Asap.compliance(::AsapOptim.OptResults, ::AsapOptim.OptParams)
Asap.axial_force(::AsapOptim.OptResults, ::AsapOptim.OptParams)
```
