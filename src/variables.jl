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
    SectionVariable <: AbstractVariable

A geometric section-property variable: the design entry REPLACES one of
the element section's properties (absolute semantics, like
[`AreaVariable`](@ref)). `property` is `:A`, `:Ix`, `:Iy`, or `:J` —
`:A` is exactly an [`AreaVariable`](@ref); the flexural/torsional
properties require a `FrameElement` (they do nothing on a truss). The
section's other properties and material are kept.

# Fields
- `element`: the controlled element (`FrameElement`, or `TrussElement` for `:A`)
- `value`, `lb`, `ub`: starting value and bounds [length⁴, or length² for `:A`]
- `property::Symbol`: `:A`, `:Ix`, `:Iy`, or `:J`

# Constructors
```julia-repl
julia> SectionVariable(beam, 8e7, 1e6, 5e8, :Ix)

julia> SectionVariable(beam, 1e6, 5e8, :Ix)   # starting value = current section value
```
"""
struct SectionVariable <: AbstractVariable
    element::Union{FrameElement,TrussElement}
    value::Float64
    lb::Float64
    ub::Float64
    property::Symbol

    function SectionVariable(element::Union{FrameElement,TrussElement},
        value::Real, lb::Real, ub::Real, property::Symbol=:A)
        property in (:A, :Ix, :Iy, :J) ||
            throw(ArgumentError("property must be :A, :Ix, :Iy, or :J"))
        property === :A || element isa FrameElement ||
            throw(ArgumentError("flexural/torsional section variables (:$property) require a FrameElement"))
        @assert element.section isa Section "section variables require a geometric Section " *
                                            "(RigiditySections have no geometric properties — vary their rigidities instead)"
        @assert 0 < lb <= value <= ub "need 0 < lb ≤ value ≤ ub"
        return new(element, value, lb, ub, property)
    end
end

SectionVariable(element::Union{FrameElement,TrussElement}, lb::Real, ub::Real,
    property::Symbol=:A) =
    SectionVariable(element, getproperty(element.section, property), lb, ub, property)

"""
    CoupledVariable <: AbstractVariable

A variable that mirrors an independent parent variable's design entry,
scaled by `factor` — e.g. symmetric node pairs (`factor = −1` for mirrored
x-coordinates), or member groups sharing one area.

By default a node coupled to a [`SpatialVariable`](@ref) moves along the
PARENT's direction: `x′ = x₀ + factor · value · d_parent`. For a spatial
parent, the coupled node can instead move along its own direction:

- `direction = d`: the target moves `factor · value` along the unit vector
  `d` (normalized at construction);
- `mirror = :X` (or `:Y`, `:Z`, or a plane-normal vector `n`): the target's
  direction is the parent's direction REFLECTED across a plane with normal
  `n`, `d = d_p − 2(d_p·n)n` — the mirror image of a symmetric pair, so one
  design entry drives both halves of a symmetric truss or shape (with
  `factor = 1`, the pair stays exactly symmetric for any value).

# Fields
- `target`: the coupled node or element
- `parent::AbstractVariable`: the independent variable whose entry is shared
- `factor::Float64`: multiplier applied to the parent's entry
- `direction::Union{Nothing,SVector{3,Float64}}`: the target's own unit
  direction (spatial couplings only), `nothing` = the parent's direction

# Examples
```julia-repl
julia> parent = SpatialVariable(left_node, 0.0, -1.0, 1.0, :X)

julia> CoupledVariable(right_node, parent, -1.0)   # mirror across the axis

julia> rail = SpatialVariable(left_panel_point, [cos(θ), sin(θ), 0.0], -0.3, 0.3)   # along a pitched chord

julia> CoupledVariable(right_panel_point, rail; mirror = :X)   # along the mirrored chord, symmetric
```
"""
struct CoupledVariable <: AbstractVariable
    target::Union{Node,FrameElement,TrussElement,Asap.FDMelement,
        Tuple{FrameElement,Symbol}}    # (element, :start/:end/:both) for joint parents
    parent::AbstractVariable
    factor::Float64
    direction::Union{Nothing,SVector{3,Float64}}

    function CoupledVariable(target, parent::AbstractVariable, factor::Real=1.0;
        direction=nothing, mirror=nothing)
        @assert !(parent isa CoupledVariable) "chain couplings to the INDEPENDENT parent variable"
        d = nothing
        if direction !== nothing || mirror !== nothing
            parent isa SpatialVariable ||
                throw(ArgumentError("direction/mirror apply only to couplings of a SpatialVariable"))
            target isa Node || throw(ArgumentError("a coupling with a direction must target a Node"))
            direction === nothing || mirror === nothing ||
                throw(ArgumentError("give either direction or mirror, not both"))
            if direction !== nothing
                d = _unit3(direction, "direction")
            else
                n = mirror isa Symbol ? _axis_vector(mirror) : _unit3(mirror, "mirror normal")
                dp = parent.direction
                d = dp - 2 * dot(dp, n) * n
            end
        end
        return new(target, parent, Float64(factor), d)
    end
end

function _unit3(v, what)
    length(v) == 3 || throw(ArgumentError("$what must be in R³"))
    u = SVector{3,Float64}(v...)
    n = norm(u)
    n > 0 || throw(ArgumentError("$what must be a nonzero vector"))
    return u / n
end

"the direction a node coupled to a spatial parent moves along (unit)"
_coupled_direction(v::CoupledVariable) = v.direction === nothing ? v.parent.direction : v.direction

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

"""
    ThicknessVariable <: AbstractVariable

A shell-thickness variable: the design entry REPLACES the thickness of every
shell element in `elements` (absolute semantics, like
[`AreaVariable`](@ref)) — one variable per wall, slab panel or group of
panels poured/cast at one thickness. Stiffness (membrane ∝ t, bending ∝ t³,
shear ∝ t) AND self-weight (∝ ρt, with design-dependent loads) follow it.

Requires `ShellSection`s (a `ShellRigiditySection` has no thickness to vary);
each element's material is kept.

# Fields
- `elements::Vector`: the controlled `ShellElement`s
- `value`, `lb`, `ub`: starting thickness and bounds [length]

# Constructors
```julia-repl
julia> ThicknessVariable(wall.elements, 0.2, 0.15, 0.4)

julia> ThicknessVariable(wall.elements, 0.15, 0.4)   # starting value = current thickness
```
"""
struct ThicknessVariable <: AbstractVariable
    elements::Vector{Any}
    value::Float64
    lb::Float64
    ub::Float64

    function ThicknessVariable(elements::AbstractVector, value::Real, lb::Real, ub::Real)
        isempty(elements) && throw(ArgumentError("a ThicknessVariable needs at least one element"))
        for el in elements
            el isa ShellElement || throw(ArgumentError("ThicknessVariable controls ShellElements only"))
            el.section isa ShellSection || throw(ArgumentError(
                "ThicknessVariable needs ShellSections (a ShellRigiditySection has no thickness)"))
        end
        @assert 0 < lb <= value <= ub "need 0 < lb ≤ value ≤ ub"
        return new(collect(Any, elements), value, lb, ub)
    end
end

ThicknessVariable(el::ShellElement, value::Real, lb::Real, ub::Real) = ThicknessVariable([el], value, lb, ub)
ThicknessVariable(elements::AbstractVector, lb::Real, ub::Real) =
    ThicknessVariable(elements, Asap.thickness(first(elements).section), lb, ub)
