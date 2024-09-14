"""
    NetworkResults

Results of a differentiable FDM analysis of a network.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- Q: force densities of all elements
"""
struct NetworkResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    Q::Vector{Float64}
end