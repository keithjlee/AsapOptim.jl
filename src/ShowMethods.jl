#=
REPL display methods for AsapOptim types (FormaSlab pattern).
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

function Base.show(io::IO, ::MIME"text/plain", v::CoupledVariable)
    tgt = v.target isa Node ? "node :$(v.target.id)" : "element :$(v.target.id)"
    println(io, "CoupledVariable  (shares its parent's design entry)")
    print(io, "  target $tgt, factor = $(v.factor)")
end

function Base.show(io::IO, ::MIME"text/plain", p::OptParams)
    println(io, "OptParams")
    println(io, "  $(length(p.values)) design variables " *
                "($(nnz(p.Sx)) position couplings, $(count(p.amask)) variable areas)")
    println(io, "  model: $(length(p.model.nodes)) nodes, $(length(p.model.elements)) elements")
    print(io, "  evaluate with solve_structure(x, p); geometry-only with GeometricProperties(x, p)")
end

function Base.show(io::IO, ::MIME"text/plain", r::OptResults)
    println(io, "OptResults")
    println(io, "  max |u| = $(maximum(abs, r.U))  [length or rad]")
    println(io, "  total volume Σ A·L = $(dot(r.A, r.L))  [length³]")
    print(io, "  fields: U (6/node), X (3×n), A, L, sections")
end
