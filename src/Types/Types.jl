# Abstract types
abstract type AbstractVariable end
abstract type IndependentVariable <: AbstractVariable end
abstract type AbstractIndexer end
abstract type AbstractOptParams end

include("Utilities.jl")

# Variables
include("Variables/Variables.jl")

# Independent structural variables
export SpatialVariable
export AreaVariable
export SectionVariable

# Independent force density variable
export QVariable

# Coupled variable
export CoupledVariable

# Unions
export TrussVariable
export NetworkVariable
export FrameVariable

# Indexers
include("Indexers/Indexers.jl")

# Optimization parameters
include("Parameters/Parameters.jl")
export TrussOptParams
export NetworkOptParams
export FrameOptParams

# results
include("Results/Results.jl")
export TrussResults
export NetworkResults
export GeometricProperties