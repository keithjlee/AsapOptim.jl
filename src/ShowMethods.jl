#=
REPL display methods for AsapOptim types (FormaSlab pattern).

Every user-facing type gets a `text/plain` show method that answers, at a
glance: what is this, what does it control/contain, and what do I do with
it next. Semantics that trip users up (ADDITIVE vs ABSOLUTE variables, the
6-DOF displacement layout) are stated right in the display.
=#

function Base.show(io::IO, ::MIME"text/plain", v::SpatialVariable)
    println(io, "SpatialVariable  (node position, ADDITIVE along one axis)")
    println(io, "  node :$(v.node.id), axis $(v.axis)")
    print(io, "  start = $(v.value), bounds [$(v.lb), $(v.ub)]  [length]")
end

function Base.show(io::IO, ::MIME"text/plain", v::AreaVariable)
    println(io, "AreaVariable  (section area, ABSOLUTE)")
    println(io, "  element :$(v.element.id)")
    print(io, "  start = $(v.value), bounds [$(v.lb), $(v.ub)]  [length²]")
end

function Base.show(io::IO, ::MIME"text/plain", v::JointVariable)
    println(io, "JointVariable  (rotational end-spring stiffness ky = kz, ABSOLUTE)")
    println(io, "  element :$(v.element.id), position $(v.position)")
    print(io, "  start = $(v.value), bounds [$(v.lb), $(v.ub)]  [force·length/rad]")
end

function Base.show(io::IO, ::MIME"text/plain", v::QVariable)
    println(io, "QVariable  (FDM force density, ABSOLUTE)")
    println(io, "  element :$(v.element.id)")
    print(io, "  start = $(v.value), bounds [$(v.lb), $(v.ub)]  [force/length]")
end

# a CoupledVariable's target can be a Node, an element (frame/truss/FDM),
# or an (element, position) tuple for joint parents — describe each legibly
_target_description(t::Node) = "node :$(t.id)"
_target_description(t::Tuple) = "element :$(t[1].id) ($(t[2]) joint)"
_target_description(t) = "element :$(t.id)"

function Base.show(io::IO, ::MIME"text/plain", v::CoupledVariable)
    println(io, "CoupledVariable  (shares its parent's design entry)")
    println(io, "  target $(_target_description(v.target)), factor = $(v.factor)")
    print(io, "  parent: $(nameof(typeof(v.parent)))")
end

_plural(n) = n == 1 ? "" : "s"

function Base.show(io::IO, ::MIME"text/plain", p::OptParams)
    njoints = count(!iszero, p.jslot1) + count(!iszero, p.jslot2)
    println(io, "OptParams")
    println(io, "  $(length(p.values)) design variable$(_plural(length(p.values))) " *
                "($(nnz(p.Sx)) position couplings, $(count(p.amask)) variable areas" *
                (njoints > 0 ? ", $njoints variable joints)" : ")"))
    println(io, "  model: $(length(p.model.nodes)) nodes, $(length(p.model.elements)) elements")
    print(io, "  evaluate with solve_structure(x, p); geometry-only with GeometricProperties(x, p)")
end

function Base.show(io::IO, ::MIME"text/plain", r::OptResults)
    println(io, "OptResults")
    println(io, "  max |u| = $(maximum(abs, r.U))  [length or rad]")
    println(io, "  total volume Σ A·L = $(dot(r.A, r.L))  [length³]")
    print(io, "  fields: U (6/node), X (3×n), A, L, sections, EA")
end

function Base.show(io::IO, ::MIME"text/plain", g::GeometricProperties)
    println(io, "GeometricProperties  (geometry only — no structural solve)")
    println(io, "  $(size(g.X, 2)) nodes, $(length(g.L)) elements")
    println(io, "  total volume Σ A·L = $(dot(g.A, g.L))  [length³]")
    print(io, "  fields: X (3×n), L, A")
end

function Base.show(io::IO, ::MIME"text/plain", p::NetworkOptParams)
    println(io, "NetworkOptParams  (force-density form-finding)")
    println(io, "  $(length(p.values)) design variable$(_plural(length(p.values))) " *
                "($(count(p.qmask)) variable force densities)")
    println(io, "  network: $(size(p.embed_free, 1)) nodes " *
                "($(size(p.embed_free, 2)) free, $(size(p.embed_fixed, 2)) fixed), " *
                "$(length(p.qmask)) elements")
    print(io, "  evaluate with solve_network(x, p)")
end

function Base.show(io::IO, ::MIME"text/plain", t::DesignTangents)
    println(io, "DesignTangents  (implicit-diff solution tangents ∂U/∂x)")
    println(io, "  dU: $(size(t.dU, 1)) dofs × $(size(t.dU, 2)) design variables")
    print(io, "  chain with axial_force_jacobian/axial_stress_jacobian, " *
              "or slice dU for displacement rows")
end

function Base.show(io::IO, ::MIME"text/plain", r::NetworkOptResults)
    println(io, "NetworkOptResults  (form-found FDM geometry)")
    println(io, "  $(length(r.X)) nodes, $(length(r.L)) elements")
    println(io, "  z ∈ [$(minimum(r.Z)), $(maximum(r.Z))], " *
                "q ∈ [$(minimum(r.Q)), $(maximum(r.Q))]")
    print(io, "  fields: X, Y, Z, Q, L; member forces via member_forces(res)")
end
