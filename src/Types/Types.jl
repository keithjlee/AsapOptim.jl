abstract type AbstractOptParams end
abstract type AbstractVariable end
abstract type TrussOptVariable end
abstract type AbstractIndexer end

# variables
include("Variables.jl")
export SpatialVariable
export AreaVariable
export CoupledVariable

# indexers
include("Indexers.jl")
export TrussOptIndexer

# parameters
include("Parameters.jl")
export TrussOptParams

# results
include("Results.jl")
export updatemodel
export TrussResults
# export OptimResults

export TrussVariable

# misc. functions
function cleartrace!(params::TrussOptParams)
    empty!(params.losstrace)
    empty!(params.valtrace)
end
export cleartrace!
