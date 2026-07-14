#=
Shared problem for the AD-backend examples: a small planar truss with mixed
geometry + sizing variables, and a compliance objective through the full
solve. Each backend file includes this, then differentiates `OBJ` its way.
=#

using Asap, AsapOptim, LinearAlgebra

function build_problem()
    steel = Material(200e6, 80.0, 0.3)
    sec = Section(steel, 1e-2)
    rot = [true, true, true]

    n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = Node([4.0, 0.0, 0.0], vcat([true, true, false], rot))
    n3 = Node([8.0, 0.0, 0.0], vcat([false, false, false], rot))
    n4 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot), :top)
    n5 = Node([6.0, 3.0, 0.0], vcat([true, true, false], rot), :top)

    els = AbstractElement{Float64}[
        TrussElement(n1, n2, sec, :chord), TrussElement(n2, n3, sec, :chord),
        TrussElement(n4, n5, sec, :chord),
        TrussElement(n1, n4, sec, :web), TrussElement(n4, n2, sec, :web),
        TrussElement(n2, n5, sec, :web), TrussElement(n5, n3, sec, :web)]
    loads = AbstractLoad{Float64}[
        NodeForce(n2, [0.0, -50.0, 0.0]), NodeForce(n5, [10.0, -30.0, 0.0])]
    model = Asap.Model([n1, n2, n3, n4, n5], els, loads)

    v_y = SpatialVariable(n4, 0.0, -1.0, 1.5, :Y)
    a_web = AreaVariable(els[4], 1e-2, 1e-4, 5e-2)
    vars = AbstractVariable[
        v_y,
        CoupledVariable(n5, v_y),                # symmetric top nodes
        a_web,
        CoupledVariable(els[5], a_web),          # grouped web areas
        CoupledVariable(els[6], a_web),
        CoupledVariable(els[7], a_web),
        AreaVariable(els[1], 1e-2, 1e-4, 5e-2)]

    params = TrussOptParams(model, vars)
    return params, copy(params.values)
end

const params, x0 = build_problem()

# objective: compliance through the full differentiable solve.
# defined as a NAMED function (not an anonymous closure over locals) — this
# matters for Enzyme, see 03_enzyme.jl
OBJ(x) = compliance(solve_structure(x, params), params)

# quick timing helper: median of n warm evaluations
med(v) = sort(v)[cld(length(v), 2)]
function time_gradient(g; n = 7)
    g()
    med([(@elapsed g()) * 1e3 for _ in 1:n])
end
