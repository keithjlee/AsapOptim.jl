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
    lb::Float64 #lower bound of variable
    ub::Float64 #upper bound of variable
    axis::Symbol #which spatial coordinate?
    iglobal::Int64 # position in the vector of active design variables

    function SpatialVariable(nodeindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        # @assert in(axis, validaxes)

        new(nodeindex, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        if typeof(node) == FDMnode
            @assert node.dof == false "FDM spatial variable only apply to anchor (fixed) nodes"
        end

        new(node.nodeID, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        # @assert in(axis, validaxes)
        if typeof(node) == FDMnode
            @assert node.dof == false "FDM spatial variable only apply to anchor (fixed) nodes"
        end

        value = node.position[axis2ind[axis]]

        new(node.nodeID, value, lowerbound, upperbound, axis)
    end
end

"""
    AreaVariable <: IndependentVariable

A variable tied to the area of an element

```julia
AreaVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
AreaVariable(element::Asap.AbstractElement, value::Float64, lowerbound::Float64, upperbound::Float64)
AreaVariable(element::Asap.AbstractElement, lowerbound::Float64, upperbound::Float64)
```
"""
mutable struct AreaVariable <: IndependentVariable
    i::Int64 #index of element, e.g. A[i] is the area variable
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Int64

    function AreaVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
        new(elementindex, value, lowerbound, upperbound)
    end

    function AreaVariable(element::Asap.AbstractElement, value::Float64, lowerbound::Float64, upperbound::Float64)
        new(element.elementID, value, lowerbound, upperbound)
    end

    function AreaVariable(element::Asap.AbstractElement, lowerbound::Float64, upperbound::Float64)
        new(element.elementID, element.section.A, lowerbound, upperbound)
    end
end

include("NetworkVariables.jl")
include("FrameVariables.jl")
include("CoupledVariables.jl")

const TrussVariable = Union{SpatialVariable, AreaVariable, CoupledVariable}
const NetworkVariable = Union{SpatialVariable, QVariable, CoupledVariable}
const FrameVariable = Union{SpatialVariable, AreaVariable, SectionVariable, CoupledVariable}
