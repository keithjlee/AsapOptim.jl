"""
    TrussOptIndexer

Translation layer between active design variables and their indices in the global structural property vectors
"""
mutable struct TrussOptIndexer <: AbstractIndexer
    iX::Vector{Int64} #index of spatial X position in vector of all X positions
    iXg::Vector{Int64} #index of associated variable in vector of design variables
    fX::Vector{Float64} #scalar factor on design variable value
    iY::Vector{Int64} #above for Y
    iYg::Vector{Int64}
    fY::Vector{Float64}
    iZ::Vector{Int64} #above for Z
    iZg::Vector{Int64}
    fZ::Vector{Float64}
    iA::Vector{Int64} #above  for Area
    iAg::Vector{Int64}
    fA::Vector{Float64}
    activeX::Bool
    activeY::Bool
    activeZ::Bool
    activeA::Bool
end


function populate!(indexer::TrussOptIndexer, var::AreaVariable)
    push!(getfield(indexer, :iA), var.i)
    push!(getfield(indexer, :iAg), var.iglobal)
    push!(getfield(indexer, :fA), 1.)

    indexer.activeA = true
end

function populate!(indexer::TrussOptIndexer, var::CoupledVariable)
    if typeof(var.referencevariable) == SpatialVariable
        field_local, field_global, field_factor = axis2field[var.referencevariable.axis]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
        push!(getfield(indexer, field_factor), var.factor)
    else
        push!(getfield(indexer, :iA), var.i)
        push!(getfield(indexer, :iAg), var.referencevariable.iglobal)
        push!(getfield(indexer, :fA), var.factor)
    end
end


"""
    TrussOptIndexer(vars::Vector{TrussVariable})

Generate the index translation layer between model parameters and design variables
"""
function TrussOptIndexer(vars::Vector{T}) where T <: TrussVariable
    indexer = TrussOptIndexer(Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        false,
        false,
        false,
        false
        )

    for var in vars
        populate!(indexer, var)
    end
    
    indexer
end


