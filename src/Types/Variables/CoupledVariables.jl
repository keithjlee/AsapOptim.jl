"""
    CoupledVariable <: AbstractVariable

A variable assigned to an element or node whose value is dependent on the values of another variable.
"""
mutable struct CoupledVariable{T<:IndependentVariable} <: AbstractVariable
    i::Int64
    target::UInt64
    factor::Float64
    iglobal::Int64
end

"""
    CoupledVariable(node::T, ref::SpatialVariable, vec::Vector{Float64}) where {T<:Asa.AbstractNode}

Make a `SpatialVariable` for `node` that is coupled to `ref.value` but with along a different vector `vector`
"""
function CoupledVariable(node::T, ref::SpatialVariable, vector::Vector{Float64}) where {T<:Asap.AbstractNode}

    @assert length(vec) == 3

    CoupledVariable{SpatialVariable}(node.nodeID, objectid(ref), vector, 0)

end

"""
    CoupledVariable(node::Asap.AbstractNode, ref::SpatialVariable, factor::Union{Float64, Vector{Float64}} = 1.0)

Make a `SpatialVariable` for `node` that is coupled to `ref.value` along the direction `factor` * `ref.vec`.

`factor` is either a scalar or a vector in R³.
"""
function CoupledVariable(node::Asap.AbstractNode, ref::SpatialVariable, factor::Union{Float64, Vector{Float64}} = 1.0)

    @assert length(factor) == 3 || length(factor) == 1

    CoupledVariable{SpatialVariable}(node.nodeID, objectid(ref), ref.vec .* factor, 0)
end

"""
    CoupledVariable(element::Asap.AbstractElement, ref::AreaVariable, factor::Float64 = 1.0)

Make an `AreaVariable` for `element` that is coupled to `factor * ref.value`

`factor` must be > 0.
"""
function CoupledVariable(element::Asap.AbstractElement, ref::AreaVariable, factor::Float64 = 1.0)

    @assert factor > 0

    CoupledVariable{AreaVariable}(element.elementID, objectid(ref), factor, 0)
end

"""
    CoupledVariable(element::Asap.Element, ref::T, factor::Float64 = 1.0) where {T<:SectionVariable}

Make an `SectionVariable` for `element` that is coupled to `factor * ref.value`

`factor` must be > 0.
"""
function CoupledVariable(element::Asap.Element, ref::T, factor::Float64 = 1.0) where {T<:SectionVariable}

    @assert factor > 0

    CoupledVariable{T}(element.elementID, objectid(ref), factor, 0)
end


"""
    CoupledVariable(element::Asap.FDMelement, ref::QVariable, factor::Float64 = 1.0)

Make an `Qvariable` for `element` that is coupled to `factor * ref.value`

`factor` must be > 0.
"""
function CoupledVariable(element::Asap.FDMelement, ref::QVariable, factor::Float64 = 1.0)

    @assert factor > 0

    CoupledVariable{QVariable}(element.elementID, objectid(ref), factor, 0)
end