"""
    AsapOptim

Structural optimization on top of Asap's differentiable core.

The v1.0 rewrite: Asap itself now provides the pure, AD-transparent
assembly/solve path (`ModelState` → `Asap.solve`), so this package no
longer re-implements analysis or hand-writes adjoints — it is a thin,
Zygote-free layer that maps a **design vector** to a **ModelState** and
back:

- design variables ([`SpatialVariable`](@ref), [`AreaVariable`](@ref),
  [`CoupledVariable`](@ref)) declare what the optimizer controls
- [`OptParams`](@ref) compiles them into sparse scatter maps
- [`solve_structure`](@ref) (aliases `solve_truss`/`solve_frame`) evaluates
  a design: positions/sections → pure solve → [`OptResults`](@ref)
- [`axial_force`](@ref)/[`axial_stress`](@ref)/[`compliance`](@ref) are
  pure functions of the results, differentiable end-to-end
- [`GeometricProperties`](@ref) gives geometry-only quantities (lengths,
  areas) without a solve
- [`updatemodel`](@ref) materializes an optimal design back into a solved
  `Asap.Model`

Bring your own AD engine and optimizer: loading Zygote (or anything using
ChainRulesCore) activates Asap's rule extension automatically.

NOTE: the force-density Network optimization path of v0.1.x is not yet
ported (`legacy_v0/` preserves the old implementation for reference).
"""
module AsapOptim

using Asap
import Asap: axial_force, compliance    # extended with OptResults methods
using LinearAlgebra
using SparseArrays
using StaticArrays

include("variables.jl")
export AbstractVariable, SpatialVariable, AreaVariable, CoupledVariable

include("parameters.jl")
export OptParams, TrussOptParams, FrameOptParams

include("results.jl")
export OptResults, solve_structure, solve_truss, solve_frame
export axial_stress                      # axial_force/compliance re-use Asap's names
export GeometricProperties

include("update.jl")
export updatemodel

include("ShowMethods.jl")

end # module AsapOptim
