"""
    AbstractVariable

Supertype of design variables. Each *independent* variable owns one entry
of the design vector `x`; a [`CoupledVariable`](@ref) shares its parent's
entry (optionally scaled), which is how symmetry and other equality
constraints are expressed without shrinking the design space by hand.
"""
abstract type AbstractVariable end

"""
    SpatialVariable <: AbstractVariable

A nodal position variable: the design entry ADDS to the node's base
coordinate along one global axis (legacy semantics preserved — a value of
`0` leaves the node at its modeled position).

# Fields
- `node::Node`: the controlled node
- `value`: starting perturbation [length]
- `lb`, `ub`: bounds on the perturbation [length]
- `axis::Symbol`: `:X`, `:Y`, or `:Z`

# Example
```julia-repl
julia> SpatialVariable(node, 0.0, -1.0, 1.0, :Z)   # node may move ±1 vertically
```
"""
struct SpatialVariable <: AbstractVariable
    node::Node
    value::Float64
    lb::Float64
    ub::Float64
    axis::Symbol

    function SpatialVariable(node::Node, value::Real, lb::Real, ub::Real, axis::Symbol)
        @assert axis in (:X, :Y, :Z) "axis must be :X, :Y, or :Z"
        @assert lb <= value <= ub "starting value must lie within [lb, ub]"
        return new(node, value, lb, ub, axis)
    end
end

"""
    AreaVariable <: AbstractVariable

A cross-section area variable: the design entry REPLACES the element
section's area (absolute semantics — the starting value is the initial
area). The section's other geometric properties (Ix, Iy, J) and material
are kept; only elements with a geometric `Section` can carry an area
variable (a `RigiditySection` has no area to vary — parameterize its
rigidities directly instead).

# Fields
- `element`: the controlled element (`FrameElement` or `TrussElement`)
- `value`: starting area [length²]
- `lb`, `ub`: bounds [length²]
"""
struct AreaVariable <: AbstractVariable
    element::Union{FrameElement,TrussElement}
    value::Float64
    lb::Float64
    ub::Float64

    function AreaVariable(element::Union{FrameElement,TrussElement},
        value::Real, lb::Real, ub::Real)
        @assert element.section isa Section "area variables require a geometric Section " *
                                            "(RigiditySections have no area — vary their rigidities instead)"
        @assert 0 < lb <= value <= ub "need 0 < lb ≤ value ≤ ub"
        return new(element, value, lb, ub)
    end
end

"""
    CoupledVariable <: AbstractVariable

A variable that mirrors an independent parent variable's design entry,
scaled by `factor` — e.g. symmetric node pairs (`factor = −1` for mirrored
x-coordinates), or member groups sharing one area.

# Fields
- `target`: the coupled node or element
- `parent::AbstractVariable`: the independent variable whose entry is shared
- `factor::Float64`: multiplier applied to the parent's entry

# Examples
```julia-repl
julia> parent = SpatialVariable(left_node, 0.0, -1.0, 1.0, :X)

julia> CoupledVariable(right_node, parent, -1.0)   # mirror across the axis
```
"""
struct CoupledVariable <: AbstractVariable
    target::Union{Node,FrameElement,TrussElement}
    parent::AbstractVariable
    factor::Float64

    function CoupledVariable(target, parent::AbstractVariable, factor::Real=1.0)
        @assert !(parent isa CoupledVariable) "chain couplings to the INDEPENDENT parent variable"
        return new(target, parent, Float64(factor))
    end
end
