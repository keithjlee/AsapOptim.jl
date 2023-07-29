"""
    TrussOptIndexer

Translation layer between active design variables and their indices in the global structural property vectors
"""
mutable struct TrussOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
end

function populate!(indexer::TrussOptIndexer, var::SpatialVariable)
    field_local, field_global = axis2field[var.axis]

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), var.iglobal)

end

function populate!(indexer::TrussOptIndexer, var::AreaVariable)
    push!(getfield(indexer, :iA), var.i)
    push!(getfield(indexer, :iAg), var.iglobal)

end

function populate!(indexer::TrussOptIndexer, var::CoupledVariable)
    if typeof(var.referencevariable) == SpatialVariable
        field_local, field_global = axis2field[var.referencevariable.axis]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
    else
        push!(getfield(indexer, :iA), var.i)
        push!(getfield(indexer, :iAg), var.referencevariable.iglobal)
    end
end

# function populate!(indexer::TrussOptIndexer, var::MirroredVariable)

#     field_local, field_global = axis2field[var.referencevariable.axis]

# end

"""
    TrussOptIndexer(vars::Vector{TrussVariable})

Generate the index translation layer between model parameters and design variables
"""
function TrussOptIndexer(vars::Vector{TrussVariable})
    indexer = TrussOptIndexer(Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}())

    for var in vars
        populate!(indexer, var)
    end
    
    indexer
end