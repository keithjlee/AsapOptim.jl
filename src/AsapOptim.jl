module AsapOptim

using Reexport

# Asap dependencies
using Asap, AsapToolkit

# Analysis dependencies
@reexport using LinearAlgebra, SparseArrays

# Optimization
@reexport using ChainRulesCore, Zygote
@reexport import Optimization 
@reexport using OptimizationNLopt
@reexport using IterativeSolvers

# Truss optimization
include("Truss/Translation.jl")

include("Truss/Types.jl")
export TrussOptProblem

# supertypes
export TrussVariable

# variables
export SpatialVariable, AreaVariable, CoupledVariable

include("Truss/Functions.jl")
export kglobal, L, Rtruss, assembleglobalK, solveU, Utruss, replacevalues, addvalues

include("Truss/Adjoints.jl")

end # module AsapOptim
