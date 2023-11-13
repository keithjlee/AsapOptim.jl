abstract type AbstractOptParams end
abstract type AbstractVariable end
abstract type TrussOptVariable end
abstract type AbstractIndexer end

include("utilities.jl")

# variables
include("Variables.jl")
export SpatialVariable
export AreaVariable
export QVariable
export CoupledVariable

# unions
export TrussVariable
export NetworkVariable

# indexers
include("Indexers.jl")
export TrussOptIndexer
export NetworkOptIndexer

# parameters
include("Parameters.jl")
export TrussOptParams
export NetworkOptParams

# results
include("Results.jl")
export updatemodel
export updatenetwork
export TrussResults
export NetworkResults
export GeometricProperties

# misc. functions
function cleartrace!(params::AbstractOptParams)
    empty!(params.losstrace)
    empty!(params.valtrace)
end
export cleartrace!