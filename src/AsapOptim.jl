module AsapOptim

using Reexport

# Asap dependencies
using Asap, AsapToolkit

# Analysis dependencies
@reexport using LinearAlgebra, SparseArrays

# Optimization
@reexport using ChainRulesCore, Zygote
@reexport import Optimization 
@reexport using OptimizationNLopt: NLopt
@reexport using IterativeSolvers

# Truss optimization
include("Truss/Translation.jl")

include("Truss/Types.jl")
export TrussOptParams

# supertypes
export TrussVariable

# variables
export SpatialVariable, AreaVariable, CoupledVariable

# output
export OptimResults

include("Truss/Functions.jl")
include("Truss/Adjoints.jl")
export kglobal, L, Rtruss, assembleglobalK, solveU, Utruss, replacevalues, addvalues

include("Truss/Utilities.jl")
export cleartrace!

include("Truss/ObjectiveFunctions.jl")
export displacement
export compliance

include("Truss/SecondaryFunctions.jl")
export Faxial



end # module AsapOptim
