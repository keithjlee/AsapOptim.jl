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

    function CoupledVariable(node::Node, ref::SpatialVariable, factor = 1.)
        new(node.nodeID, ref, factor)
    end

    function CoupledVariable(element::Union{Element, TrussElement}, ref::AreaVariable, factor = 1.)
        @assert factor > 0 "Coupling factor must be greater than 0 when referring to area variables"
        new(element.elementID, ref, factor)
    end

    function CoupledVariable(element::Element, ref::SectionVariable, factor = 1.)
        @assert factor > 0 "Coupling factor must be greater than 0 when referring to area variables"
        new(element.elementID, ref, factor)
    end

    function CoupledVariable(element::FDMelement, ref::QVariable, factor = 1.)
        @assert factor > 0 "Coupling factor must be greater than 0 when referring to force density variables"
        new(element.elementID, ref, factor)
    end
end