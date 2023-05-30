abstract type FrameOptVariable end

mutable struct FrameOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
    iIz::Vector{Int64}
    iIzg::Vector{Int64}
    iIy::Vector{Int64}
    iIyg::Vector{Int64}
end