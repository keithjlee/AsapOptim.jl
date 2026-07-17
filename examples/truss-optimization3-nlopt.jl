#=
Example 5 — constrained minimum volume, NLopt driven DIRECTLY, with the
AD engine chosen PER FUNCTION:

  * objective (material volume): reverse mode (Zygote) — one scalar output
  * constraints (25 displacement + 47 stress rows): FORWARD mode
    (ForwardDiff, prepared) — a full-Jacobian workload where forward AD
    beats reverse by two orders of magnitude (see docs/AD_BACKENDS.md)
  * one Asap.CachedSolver on the params: the objective gradient, the
    constraint Jacobian, and every ForwardDiff chunk at a design iterate
    share a single stiffness factorization

Same Warren truss and constraint set as example 2 (truss-optimization2.jl),
so the two are directly comparable — this one talks to NLopt's C API
conventions with no bridge packages in between.

Additional packages:  ] add NLopt Zygote ForwardDiff DifferentiationInterface
=#

using AsapOptim, Asap, LinearAlgebra
using NLopt, Zygote, ForwardDiff
using DifferentiationInterface

const VISUALIZE = get(ENV, "ASAP_EXAMPLE_VIZ", "true") == "true"
VISUALIZE && using CairoMakie

# section and constants (identical to example 2)
begin
    steel = Material(200e6, 80.0, 0.3)
    fy = 350e3                                # yield stress [kN/m²]
    section = Section(steel, 0.01)

    L = 10.0
    nx = 12
    dy = 1.0
    load = 75.0
    dx = L / nx
    dmax = L / 360
end

# generate the truss
begin
    bottom_nodes = [Node([x, 0.0, 0.0], :free, :bottom) for x in range(0, L, nx + 1)]
    fixnode!(first(bottom_nodes), :pinned)
    first(bottom_nodes).id = :pin
    fixnode!(last(bottom_nodes), :xfree)
    last(bottom_nodes).id = :roller

    top_nodes = [Node([dx / 2 + dx * (i - 1), dy, 0.0], :free, :top) for i in 1:nx]

    bottom_elements = [TrussElement(bottom_nodes[i], bottom_nodes[i+1], section, :bottom) for i in 1:nx]
    top_elements = [TrussElement(top_nodes[i], top_nodes[i+1], section, :top) for i in 1:nx-1]
    web_elements = [
        [TrussElement(nb, nt, section, :web) for (nb, nt) in zip(bottom_nodes[1:end-1], top_nodes)];
        [TrussElement(nt, nb, section, :web) for (nt, nb) in zip(top_nodes, bottom_nodes[2:end])]
    ]

    nodes = [bottom_nodes; top_nodes]
    elements = AbstractElement{Float64}[bottom_elements; top_elements; web_elements]
    loads = AbstractLoad{Float64}[NodeForce(node, [0.0, -load, 0.0]) for node in nodes[:bottom]]

    model = Asap.Model(nodes, elements, loads)
    planarize!(model)
    solve!(model)
end

element_indices = vcat([[e.nodeStart.index, e.nodeEnd.index] for e in model.elements]...)

#=
Variables and params — note the CachedSolver
=#
begin
    xmin, xmax = (-1, 1) .* dx ./ 2 .* 0.9
    ymin, ymax = -0.9dy, dy

    vars = AbstractVariable[
        [SpatialVariable(node, 0.0, xmin, xmax, :X) for node in model.nodes[:top]];
        [SpatialVariable(node, 0.0, ymin, ymax, :Y) for node in model.nodes[:top]];
        [AreaVariable(element, element.section.A, 0.01 * element.section.A, 5 * element.section.A)
         for element in model.elements]
    ]

    params = OptParams(model, vars; solver = Asap.CachedSolver())
    x0 = copy(params.values)
end

#=
Objective: material volume — geometry only, reverse-mode gradient.
NLopt convention: f(x, grad) returns the value and fills grad IN PLACE
when NLopt asks for it (length(grad) > 0).
=#
volume(x) = begin
    geo = GeometricProperties(x, params)
    dot(geo.L, geo.A)
end

function nlopt_objective(x::Vector, grad::Vector)
    if length(grad) > 0
        v, g = Zygote.withgradient(volume, x)
        grad .= g[1]
        return v
    end
    return volume(x)
end

#=
Constraints: displacement + stress rows, each normalized by its limit
(quantity/limit − 1 ≤ 0 — mixed-scale rows stall MMA, see example 2).
The Jacobian is FORWARD-mode with a prepared DifferentiationInterface
operator; NLopt's vector-constraint callback receives grad as an
(n_vars × n_constraints) matrix — the Jacobian TRANSPOSED.
=#
function constraints(x)
    res = solve_truss(x, params)
    return [
        abs.(res.U[2:6:end]) ./ dmax .- 1.0;
        abs.(axial_stress(res, params)) ./ fy .- 1.0
    ]
end

const CSTR_BACKEND = AutoForwardDiff()
const CSTR_PREP = prepare_jacobian(constraints, CSTR_BACKEND, x0)

function nlopt_constraints!(result::Vector, x::Vector, grad::Matrix)
    if length(grad) > 0
        v, J = value_and_jacobian(constraints, CSTR_PREP, CSTR_BACKEND, x)
        result .= v
        grad .= J'                       # NLopt wants (n, m)
    else
        result .= constraints(x)
    end
    return
end

n_constraints = length(constraints(x0))
@assert all(constraints(x0) .<= 0) "start from a feasible design"

#=
Optimize — NLopt MMA, no bridge packages
=#
opt = Opt(:LD_MMA, length(x0))
NLopt.lower_bounds!(opt, params.lb)
NLopt.upper_bounds!(opt, params.ub)
NLopt.min_objective!(opt, nlopt_objective)
NLopt.inequality_constraint!(opt, nlopt_constraints!, fill(1e-8, n_constraints))
NLopt.maxeval!(opt, 500)
NLopt.maxtime!(opt, 60)

v0 = volume(x0)
stats = @timed optimize(opt, x0)
minf, minx, ret = stats.value
c_final = constraints(minx)

println("initial volume:   ", v0, " m³")
println("optimized volume: ", minf, " m³  (", round(minf / v0 * 100; digits=1), "% of initial)")
println("NLopt status: ", ret, " after ", NLopt.numevals(opt), " evaluations (",
    round(stats.time; digits=2), " s)")
println("constraints satisfied: ", all(c_final .<= 1e-6),
    "  (max = ", maximum(c_final), ")")
solver = params.solver
println("factorizations shared: ", solver.hits, " cache hits, ",
    solver.refactorizations, " numeric refactorizations")

model2 = updatemodel(params, minx)

if VISUALIZE
    p2 = Point2.([n.position[1:2] for n in model2.nodes])
    e2 = p2[element_indices]

    areas = [e.section.A for e in model2.elements]
    lw = areas ./ maximum(areas) .* 4

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    ylims!(-dy / 2, 3dy)
    hidespines!(ax); hidedecorations!(ax)
    linesegments!(e2, color = :black, linewidth = lw)
    scatter!(p2, color = :white, strokecolor = :black, strokewidth = 1)
    display(fig)
end
