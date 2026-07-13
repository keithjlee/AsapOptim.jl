include("TrussIndexer.jl")
include("NetworkIndexer.jl")
include("FrameIndexer.jl")

function populate!(indexer::AbstractIndexer, var::SpatialVariable)

    x, y, z = var.vec

    if x != 0.0
        push!(indexer.iX, var.i)
        push!(indexer.iXg, var.iglobal)
        push!(indexer.fX, x)

        indexer.activeX = true
    end

    if y != 0.0
        push!(indexer.iY, var.i)
        push!(indexer.iYg, var.iglobal)
        push!(indexer.fY, y)

        indexer.activeY = true
    end

    if z != 0.0
        push!(indexer.iZ, var.i)
        push!(indexer.iZg, var.iglobal)
        push!(indexer.fZ, z)

        indexer.activeZ = true
    end

end

function populate!(indexer::Union{TrussOptIndexer, FrameOptIndexer}, var::AreaVariable)

    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, 1.0)

    indexer.activeA = true

end

function populate!(indexer::AbstractIndexer, var::CoupledVariable{SpatialVariable})

    fx, fy, fz = var.factor

    if fx != 0.0
        push!(indexer.iX, var.i)
        push!(indexer.iXg, var.iglobal)
        push!(indexer.fX, fx)

        indexer.activeX = true
    end

    if fy != 0.0
        push!(indexer.iY, var.i)
        push!(indexer.iYg, var.iglobal)
        push!(indexer.fY, fy)

        indexer.activeY = true
    end

    if fz != 0.0
        push!(indexer.iZ, var.i)
        push!(indexer.iZg, var.iglobal)
        push!(indexer.fZ, fz)

        indexer.activeZ = true
    end
end

function populate!(indexer::Union{TrussOptIndexer, FrameOptIndexer}, var::CoupledVariable{AreaVariable})

    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionA}})

    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionIx}})

    push!(indexer.iIx, var.i)
    push!(indexer.iIxg, var.iglobal)
    push!(indexer.fIx, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionIy}})

    push!(indexer.iIy, var.i)
    push!(indexer.iIyg, var.iglobal)
    push!(indexer.fIy, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionJ}})

    push!(indexer.iJ, var.i)
    push!(indexer.iJg, var.iglobal)
    push!(indexer.fJ, var.factor)
    
end