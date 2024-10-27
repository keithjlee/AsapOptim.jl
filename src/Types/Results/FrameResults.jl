"""
    TrussResults

Results of a differentiable structural analysis of a truss model.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- A: Area of all elements
- L: Length of all elements
- K: Elemental stiffness matrices in GCS
- R: Elemental transformation matrices
- U: Displacement vector of all nodes
"""
struct FrameResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    A::Vector{Float64}
    Ix::Vector{Float64}
    Iy::Vector{Float64}
    J::Vector{Float64}
    L::Vector{Float64}
    K::Vector{Matrix{Float64}}
    R::Vector{Matrix{Float64}}
    U::Vector{Float64}
end

