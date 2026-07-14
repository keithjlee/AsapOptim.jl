#=
Example 1 — geometry optimization of a Warren-style truss (compliance)

The classic AsapOptim demo, on the v1.0 API: the top-chord nodes of a
planar truss are free to move in x and y; find the geometry of minimum
compliance (maximum stiffness) under gravity loads.

Additional packages:  ] add Nonconvex NonconvexNLopt Zygote CairoMakie
(kjlMakie optional, for Keith's plot theme)
=#

using AsapOptim, Asap, LinearAlgebra
using Nonconvex, NonconvexNLopt, Zygote

# set ASAP_EXAMPLE_VIZ=false to run headless (e.g. for testing)
const VISUALIZE = get(ENV, "ASAP_EXAMPLE_VIZ", "true") == "true"
VISUALIZE && using CairoMakie

# section and problem constants
begin
    steel = Material(200e6, 80.0, 0.3)        # E [kN/m²], ρ, ν
    section = Section(steel, 0.01)            # axial-only, A = 0.01 m²

    L = 10.0    # total span
    nx = 12     # bays
    dy = 1.0    # truss depth
    load = 75.0 # per-node load [kN]
end

# generate the truss
begin
    dx = L / nx

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

# node pairs for plotting element segments
element_indices = vcat([[e.nodeStart.index, e.nodeEnd.index] for e in model.elements]...)

if VISUALIZE
    p = Point2.([n.position[1:2] for n in model.nodes])
    e = p[element_indices]

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    ylims!(-dy / 2, 2dy)
    hidespines!(ax); hidedecorations!(ax)
    linesegments!(e, color = :black)
    scatter!(p, color = :white, strokecolor = :black, strokewidth = 1)
    display(fig)
end

#=
Spatial optimization — compliance
=#

# spatial variables for all top nodes (ADDITIVE perturbations)
xmin, xmax = (-1, 1) .* dx ./ 2 .* 0.9
ymin, ymax = -0.9dy, dy

vars = AbstractVariable[
    [SpatialVariable(node, 0.0, xmin, xmax, :X) for node in model.nodes[:top]];
    [SpatialVariable(node, 0.0, ymin, ymax, :Y) for node in model.nodes[:top]]
]

params = TrussOptParams(model, vars)
x0 = copy(params.values)

# objective: compliance Fᵀu (external work — lower is stiffer)
OBJ = x -> compliance(solve_truss(x, params), params)

# test value and gradient (gradients flow through Asap's pure solve)
o0, ∇o0 = Zygote.withgradient(OBJ, x0)
println("initial compliance: ", o0)

# optimize with L-BFGS via Nonconvex
optmodel = Nonconvex.Model(OBJ)
addvar!(optmodel, params.lb, params.ub)

alg = NLoptAlg(:LD_LBFGS)
opts = NLoptOptions(maxeval = 500, maxtime = 60)

res = optimize(optmodel, alg, x0, options = opts)
println("optimized compliance: ", res.minimum, "  (", round(res.minimum / o0 * 100; digits=1), "% of initial)")

# materialize the optimal design as a solved Asap model
model2 = updatemodel(params, res.minimizer)

if VISUALIZE
    p2 = Point2.([n.position[1:2] for n in model2.nodes])
    e2 = p2[element_indices]

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    ylims!(-dy / 2, 3dy)
    hidespines!(ax); hidedecorations!(ax)
    linesegments!(e2, color = :black)
    scatter!(p2, color = :white, strokecolor = :black, strokewidth = 1)
    display(fig)
end
