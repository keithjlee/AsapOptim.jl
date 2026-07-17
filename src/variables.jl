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
position along a direction — a global axis, or an arbitrary directional
"rail" (legacy semantics preserved — a value of `0` leaves the node at
its modeled position):

    x′ = x₀ + value · direction

The rail direction is NORMALIZED at construction, so `value`/`lb`/`ub`
are arc-length along the rail [length] regardless of the vector you pass.

# Fields
- `node::Node`: the controlled node
- `value`: starting perturbation along the direction [length]
- `lb`, `ub`: bounds on the perturbation [length]
- `direction::SVector{3,Float64}`: unit direction of travel

# Constructors
```julia-repl
julia> SpatialVariable(node, 0.0, -1.0, 1.0, :Z)          # global axis (:X/:Y/:Z)

julia> SpatialVariable(node, [1.0, 1.0, 0.0], 0.0, -1.0, 1.0)  # rail: ±1 along the diagonal

julia> SpatialVariable(node, [1.0, 1.0, 0.0], -1.0, 1.0)       # rail, starting value 0
```
"""
struct SpatialVariable <: AbstractVariable
    node::Node
    value::Float64
    lb::Float64
    ub::Float64
    direction::SVector{3,Float64}

    function SpatialVariable(node::Node, value::Real, lb::Real, ub::Real,
        direction::SVector{3,Float64})
        @assert lb <= value <= ub "starting value must lie within [lb, ub]"
        n = norm(direction)
        @assert n > 0 "direction must be a nonzero vector"
        return new(node, value, lb, ub, direction / n)
    end
end

const _AXIS_VECTORS = (X = SVector(1.0, 0.0, 0.0),
    Y = SVector(0.0, 1.0, 0.0), Z = SVector(0.0, 0.0, 1.0))

function _axis_vector(axis::Symbol)
    a = Symbol(uppercase(String(axis)))
    haskey(_AXIS_VECTORS, a) || throw(ArgumentError("axis must be :X, :Y, or :Z"))
    return _AXIS_VECTORS[a]
end

SpatialVariable(node::Node, value::Real, lb::Real, ub::Real, axis::Symbol) =
    SpatialVariable(node, value, lb, ub, _axis_vector(axis))

SpatialVariable(node::Node, direction::AbstractVector{<:Real}, value::Real, lb::Real, ub::Real) =
    (length(direction) == 3 || throw(ArgumentError("direction must be in R³"));
    SpatialVariable(node, value, lb, ub, SVector{3,Float64}(direction)))

SpatialVariable(node::Node, direction::AbstractVector{<:Real}, lb::Real, ub::Real) =
    SpatialVariable(node, direction, 0.0, lb, ub)

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
    target::Union{Node,FrameElement,TrussElement,Asap.FDMelement,
        Tuple{FrameElement,Symbol}}    # (element, :start/:end/:both) for joint parents
    parent::AbstractVariable
    factor::Float64

    function CoupledVariable(target, parent::AbstractVariable, factor::Real=1.0)
        @assert !(parent isa CoupledVariable) "chain couplings to the INDEPENDENT parent variable"
        return new(target, parent, Float64(factor))
    end
end

"""
    JointVariable <: AbstractVariable

A semi-rigid connection stiffness variable: the design entry REPLACES the
rotational end-spring stiffness (`ky` = `kz`, both bending planes) of a
frame element's connection at one or both ends [force·length/rad].

Stiffer joints cost more to fabricate (welding, bolts, embedments) — this
variable is what lets an optimizer trade connection cost against structural
performance. Axial and torsional connection stiffnesses keep the element's
existing values.

# Fields
- `element::FrameElement`: the element whose connection is controlled
- `position::Symbol`: `:start`, `:end`, or `:both`
- `value`, `lb`, `ub`: starting stiffness and bounds [force·length/rad]

# Example
```julia-repl
julia> JointVariable(beam, :both, 1e5, 1e3, 1e8)
```
"""
struct JointVariable <: AbstractVariable
    element::FrameElement
    position::Symbol
    value::Float64
    lb::Float64
    ub::Float64

    function JointVariable(element::FrameElement, position::Symbol,
        value::Real, lb::Real, ub::Real)
        @assert position in (:start, :end, :both) "position must be :start, :end, or :both"
        @assert 0 < lb <= value <= ub "need 0 < lb ≤ value ≤ ub"
        return new(element, position, value, lb, ub)
    end
end
