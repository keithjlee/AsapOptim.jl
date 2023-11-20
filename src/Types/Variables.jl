"""
    SpatialVariable <: AbstractVariable

A variable tied to the spatial position of a node in a single axis.

```julia
SpatialVariable(nodeindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
```
"""
mutable struct SpatialVariable <: AbstractVariable
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
    AreaVariable <: AbstractVariable

A variable tied to the area of an element

```julia
AreaVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
AreaVariable(element::Asap.AbstractElement, value::Float64, lowerbound::Float64, upperbound::Float64)
AreaVariable(element::Asap.AbstractElement, lowerbound::Float64, upperbound::Float64)
```
"""
mutable struct AreaVariable <: AbstractVariable
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

"""
    Qvariable <: AbstractVariable

A variable tied to the force density of an FDM element

```julia
QVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
QVariable(element::Asap.AbstractElement, value::Float64, lowerbound::Float64, upperbound::Float64)
QVariable(element::Asap.AbstractElement, lowerbound::Float64, upperbound::Float64)
```
"""
mutable struct QVariable <: AbstractVariable
    i::Int64
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Int64

    function QVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
        new(elementindex, value, lowerbound, upperbound)
    end

    function QVariable(element::FDMelement, value::Float64, lowerbound::Float64, upperbound::Float64)
        new(element.elementID, value, lowerbound, upperbound)
    end

    function QVariable(element::FDMelement, lowerbound::Float64, upperbound::Float64)
        new(element.elementID, value, lowerbound, upperbound)
    end
end

"""
    CoupledVariable <: AbstractVariable

A variable whose value points to an existing variable.

```julia
var = CoupledVariable(variable::Union{TrussNode, TrussElement, FDMelement}, reference::Union{SpatialVariable, AreaVariable, QVaraible}, factor = 1.0)
```

Where `factor` is a scalar factor applied to the value of `reference` before assigning to `variable`. IE to enforce a mirrored value, use `factor = -1.0`.
"""
mutable struct CoupledVariable <: AbstractVariable
    i::Int64
    referencevariable::AbstractVariable
    factor::Real

    function CoupledVariable(node::TrussNode, ref::SpatialVariable, factor = 1.)
        new(node.nodeID, ref, factor)
    end

    function CoupledVariable(element::TrussElement, ref::AreaVariable, factor = 1.)
        @assert factor > 0 "Coupling factor must be greater than 0 when referring to area variables"
        new(element.elementID, ref, factor)
    end

    function CoupledVariable(element::FDMelement, ref::QVariable, factor = 1.)
        @assert factor > 0 "Coupling factor must be greater than 0 when referring to force density variables"
        new(element.elementID, ref, factor)
    end
end


const TrussVariable = Union{SpatialVariable, AreaVariable, CoupledVariable}
const NetworkVariable = Union{SpatialVariable, QVariable, CoupledVariable}