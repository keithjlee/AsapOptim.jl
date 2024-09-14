# convert an axis symbol to the canonical directional vectors
const axis_to_vector = Dict(
    :x => Asap.globalX,
    :X => Asap.globalX,
    :y => Asap.globalY,
    :Y => Asap.globalY,
    :z => Asap.globalZ,
    :Z => Asap.globalZ
)

"""
    SpatialVariable <: IndependentVariable

A variable tied to the spatial position of a node in a single axis.

```julia
SpatialVariable(nodeindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
```
"""
mutable struct SpatialVariable <: IndependentVariable
    i::Int64 #index of node, e.g. X[i] is the spatial variable
    val::Float64 #value
    vec::Vector{Float64} #directional vector
    lb::Float64 #lower bound of variable
    ub::Float64 #upper bound of variable
    iglobal::Int64 # position in the vector of active design variables
end

"""
    SpatialVariable(node::T, vector::Vector{Float64}, value::Float64, lowerbound::Float64, upperbound::Float64) where {T<:Asap.AbstractNode}

Define a `SpatialVariable` for a given node, in the normalized direction of `vector`, such that:
    x' = x_0 + `value` * `vector[1]`  
    y' = y_0 + `value` * `vector[2]`  
    z' = z_0 + `value` + `vector[3]`

Where: `lowerbound` ≤ `value` ≤ `upperbound`
"""
function SpatialVariable(node::T, vector::Vector{Float64}, value::Float64, lowerbound::Float64, upperbound::Float64) where {T<:Asap.AbstractNode}
    length(vector) != 3 && error("Direction vector must be in R³")
    return SpatialVariable(node.nodeID, value, normalize(vector), lowerbound, upperbound, 0)
end

"""
    SpatialVariable(node::T, vector::Vector{Float64}, lowerbound::Float64, upperbound::Float64) where {T<:Asap.AbstractNode}

Define a `SpatialVariable` for a given node, in the normalized direction of `vector`, such that:
    x' = x_0 + `value` * `vector[1]`  
    y' = y_0 + `value` * `vector[2]`  
    z' = z_0 + `value` + `vector[3]`

Where: `lowerbound` ≤ `value` ≤ `upperbound`

And `value` is set to 0.
"""
function SpatialVariable(node::T, vector::Vector{Float64}, lowerbound::Float64, upperbound::Float64) where {T<:Asap.AbstractNode}
    return SpatialVariable(node, normalize(vector), 0.0, lowerbound, upperbound)
end

"""
    SpatialVariable(node::T, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z) where {T<:Asap.AbstractNode}

Define a `SpatialVariable` for a given node, in a given Cartesian direction, `axis` ∈ [:X, :Y, :Z]
"""
function SpatialVariable(node::T, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z) where {T<:Asap.AbstractNode}
    !in(axis, keys(axis_to_vector)) && error("Axis value not recognized, choose from (:X, :Y, :Z) or define an explict vector in R³")

    vector = axis_to_vector[axis]
    return SpatialVariable(node, vector, value, lowerbound, upperbound)
end

"""
    SpatialVariable(node::T, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z) where {T<:Asap.AbstractNode}

Define a `SpatialVariable` for a given node, in a given Cartesian direction, `axis` ∈ [:X, :Y, :Z] and default value is set to 0.
"""
function SpatialVariable(node::T, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z) where {T<:Asap.AbstractNode}

    return SpatialVariable(node, 0.0, lowerbound, upperbound, axis)
end