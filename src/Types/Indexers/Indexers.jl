function populate!(indexer::AbstractIndexer, var::SpatialVariable)
    field_local, field_global, field_factor = axis2field[var.axis]

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), var.iglobal)
    push!(getfield(indexer, field_factor), 1.)
end

# quick reference to relevant field in TrussOptIndexer from variables
const axis2field = Dict(:X => (:iX, :iXg, :fX),
    :x => (:iX, :iXg, :fX),
    :Y => (:iY, :iYg, :fY),
    :y => (:iY, :iYg, :fY),
    :Z => (:iZ, :iZg, :fZ),
    :z => (:iZ, :iZg, :fZ))

const property2field = Dict(
    :A => (:iA, iAg, :fA),
    :Ix => (:iIx, iIxg, :fIx),
    :Iy => (:iIy, iIyg, :fIy),
    :J => (:iJ, iJg, :fJ)
)

include("TrussIndexer.jl")
include("NetworkIndexer.jl")
include("FrameIndexer.jl")