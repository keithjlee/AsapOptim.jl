"""
    CoupledVariable <: AbstractVariable

A variable whose value points to an existing variable.

```julia
var = CoupledVariable(variable::Union{TrussNode, TrussElement, FDMelement}, reference::Union{SpatialVariable, AreaVariable, QVaraible}, factor = 1.0)
```

Where `factor` is a scalar factor applied to the value of `reference` before assigning to `variable`. IE to enforce a mirrored value, use `factor = -1.0`.
"""
mutable struct CoupledVariable{T<:IndependentVariable} <: AbstractVariable
    i::Int64
    target::UInt64
    factor::Real
    iglobal::Int64
end

function CoupledVariable(node::Asap.AbstractNode, ref::SpatialVariable, factor::Float64 = 1.0)

    T = typeof(ref)
    return CoupledVariable{T}(node.nodeID, objectid(ref), factor, 0)
end

function CoupledVariable(element::Asap.AbstractElement, ref::AreaVariable, factor::Float64 = 1.0)

    @assert factor > 0

    T = typeof(ref)
    return CoupledVariable{T}(element.elementID, objectid(ref), factor, 0)
end

function CoupledVariable(element::Asap.Element, ref::SectionVariable, factor::Float64 = 1.0)

    @assert factor > 0

    T = typeof(ref)
    return CoupledVariable{T}(element.elementID, objectid(rf), factor, 0)
end

function CoupledVariable(element::Asap.FDMelement, ref::QVariable, factor::Float64 = 1.0)

    @assert factor > 0

    T = typeof(ref)
    return CoupledVariable{T}(element.elementID, objectid(ref), factor, 0)
end