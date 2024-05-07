mutable struct FrameOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    fX::Vector{Float64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    fY::Vector{Float64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    fZ::Vector{Float64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
    fA::Vector{Float64}
    iIx::Vector{Int64}
    iIxg::Vector{Int64}
    fIx::Vector{Float64}
    iIy::Vector{Int64}
    iIyg::Vector{Int64}
    fIy::Vector{Float64}
    iJ::Vector{Int64}
    iJg::Vector{Int64}
    fJ::Vector{Float64}
    activeX::Bool
    activeY::Bool
    activeZ::Bool
    activeA::Bool
    activeIx::Bool
    activeIy::Bool
    activeJ::Bool
end

function populate!(indexer::FrameOptIndexer, var::AreaVariable)
    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, 1.0)

    indexer.activeA = true
end

function populate!(indexer::FrameOptIndexer, var::SectionVariable)
    field_local, field_global, field_factor = property2field[var.property]

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), var.iglobal)
    push!(getfield(indexer, field_factor), 1.)

    setfield!(indexer, property2active[var.property], true)
end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable)

    #for spatial variables, find proper axis of reference
    if typeof(var.referencevariable) == SpatialVariable

        field_local, field_global, field_factor = axis2field[var.referencevariable.axis]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
        push!(getfield(indexer, field_factor), var.factor)

    elseif typeof(var.referencevariable) == SectionVariable

        field_local, field_global, field_factor = property2field[var.referencevariable.property]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
        push!(getfield(indexer, field_factor), var.factor)

    else

        push!(getfield(indexer, :iA), var.i)
        push!(getfield(indexer, :iAg), var.referencevariable.iglobal)
        push!(getfield(indexer, :fA), var.factor)
    end
end

function FrameOptIndexer(vars::Vector{T}) where T<:FrameVariable
    indexer = FrameOptIndexer(
        Vector{Int64}(),
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
        false,
        false,
        false,
        false
        )

    # iX = Vector{Int64}()
    # iXg = Vector{Int64}()
    # fX = Vector{Float64}()
    # iY = Vector{Int64}()
    # iYg = Vector{Int64}()
    # fY = Vector{Float64}()
    # iZ = Vector{Int64}()
    # iZg = Vector{Int64}()
    # fZ = Vector{Float64}()
    # iA = Vector{Int64}()
    # iAg = Vector{Int64}()
    # fA = Vector{Float64}()
    # iIx = Vector{Int64}()
    # iIxg = Vector{Int64}()
    # fIx = Vector{Float64}()
    # iIy = Vector{Int64}()
    # iIyg = Vector{Int64}()
    # fIy = Vector{Float64}()
    # iJ = Vector{Int64}()
    # iJg = Vector{Int64}()
    # fJ = Vector{Float64}()
    # activeX = false
    # activeY = false
    # activeZ = false
    # activeA = false
    # activeIx = false
    # activeIy = false
    # activeJ = false

    for var in vars
        populate!(indexer, var)
    end

    indexer
end