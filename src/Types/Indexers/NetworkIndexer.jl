mutable struct NetworkOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    fX::Vector{<:Real}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    fY::Vector{<:Real}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    fZ::Vector{<:Real}
    iQ::Vector{Int64}
    iQg::Vector{Int64}
    fQ::Vector{<:Real}
    iN::Vector{<:Real}
    activeX::Bool
    activeY::Bool
    activeZ::Bool
    activeQ::Bool
end

function populate!(indexer::NetworkOptIndexer, var::QVariable)
    push!(getfield(indexer, :iQ), var.i)
    push!(getfield(indexer, :iQg), var.iglobal)
    push!(getfield(indexer, :fQ), 1.)

    indexer.activeQ = true
end

function populate!(indexer::NetworkOptIndexer, var::CoupledVariable)
    if typeof(var.referencevariable) == SpatialVariable
        field_local, field_global, field_factor = axis2field[var.referencevariable.axis]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
        push!(getfield(indexer, field_factor), var.factor)
    else
        push!(getfield(indexer, :iQ), var.i)
        push!(getfield(indexer, :iQg), var.referencevariable.iglobal)
        push!(getfield(indexer, :fQ), var.factor)
    end
end

"""
    NetworkOptIndexer(vars::Vector{NetworkVariable})

Generate the index translation layer between network parameters and design variables
"""
function NetworkOptIndexer(vars::Vector{NetworkVariable})
    indexer = NetworkOptIndexer(
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Real}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Real}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Real}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Real}(),
        Vector{Real}(),
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