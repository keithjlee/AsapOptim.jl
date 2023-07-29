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

    function SpatialVariable(node::Asap.AbstractNode, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        # @assert in(axis, validaxes)

        new(node.nodeID, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::Asap.AbstractNode, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        # @assert in(axis, validaxes)

        value = node.position[axis2ind[axis]]

        new(node.nodeID, value, lowerbound, upperbound, axis)
    end
end

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

mutable struct CoupledVariable <: AbstractVariable
    i::Int64
    referencevariable::AbstractVariable

    function CoupledVariable(node::TrussNode, ref::SpatialVariable)
        new(node.nodeID, ref)
    end

    function CoupledVariable(element::TrussElement, ref::AreaVariable)
        new(element.elementID, ref)
    end
end

mutable struct MirroredVariable <: AbstractVariable
    i::Int64
    referencevariable::SpatialVariable

    function MirroredVariable(node::TrussNode, ref::SpatialVariable)
        new(node.nodeID, ref)
    end
end


const TrussVariable = Union{SpatialVariable, AreaVariable, CoupledVariable, MirroredVariable}