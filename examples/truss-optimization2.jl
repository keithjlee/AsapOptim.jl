#=
Example 2 — constrained minimum-volume truss (geometry + sizing)

Same Warren-style truss as example 1, now with node positions AND member
areas as design variables: minimize material volume Σ AᵢLᵢ subject to
displacement and stress constraints (the publication S4.1 problem type).

Additional packages:  ] add Nonconvex NonconvexNLopt Zygote CairoMakie
=#

using AsapOptim, Asap, LinearAlgebra
using Nonconvex, NonconvexNLopt, Zygote

const VISUALIZE = get(ENV, "ASAP_EXAMPLE_VIZ", "true") == "true"
VISUALIZE && using CairoMakie

# section and constants
begin
    steel = Material(200e6, 80.0, 0.3)
    fy = 350e3                                # yield stress [kN/m²]
    section = Section(steel, 0.01)

    L = 10.0
    nx = 12
    dy = 1.0
    load = 75.0
end

# generate the truss (identical to example 1)
begin
    dx = L / nx
    dmax = L / 360                            # displacement limit

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
Constrained minimum volume
=#

xmin, xmax = (-1, 1) .* dx ./ 2 .* 0.9
ymin, ymax = -0.9dy, dy

vars = AbstractVariable[
    [SpatialVariable(node, 0.0, xmin, xmax, :X) for node in model.nodes[:top]];
    [SpatialVariable(node, 0.0, ymin, ymax, :Y) for node in model.nodes[:top]];
    [AreaVariable(element, element.section.A, 0.01 * element.section.A, 5 * element.section.A)
     for element in model.elements]
]

params = TrussOptParams(model, vars)
x0 = copy(params.values)

# objective: material volume — geometry-only, no solve needed
function obj(x, p)
    geo = GeometricProperties(x, p)
    dot(geo.L, geo.A)
end
OBJ = x -> obj(x, params)
o0, ∇o0 = Zygote.withgradient(OBJ, x0)
println("initial volume: ", o0, " m³")

# constraints: vertical displacements and axial stresses, each row
# NORMALIZED by its limit (quantity/limit − 1 ≤ 0). Normalization matters:
# raw displacement rows are O(0.01) while stress rows are O(1e5), and that
# scale mismatch stalls MMA at an infeasible point with ~2× the volume.
# (NOTE the v1.0 DOF layout: 6 DOFs per node — vertical is U[2:6:end])
function cstr(x, p, dmax, smax)
    res = solve_truss(x, p)

    vertical_displacements = res.U[2:6:end]
    stresses = axial_stress(res, p)

    return [
        abs.(vertical_displacements) ./ dmax .- 1.0;
        abs.(stresses) ./ smax .- 1.0
    ]
end
CSTR = x -> cstr(x, params, dmax, fy)
c0, ∇c0 = Zygote.withjacobian(CSTR, x0)
@assert all(c0 .<= 0) "start from a feasible design (increase the initial area if this fires)"

# optimize with MMA
optmodel = Nonconvex.Model(OBJ)
addvar!(optmodel, params.lb, params.ub)
add_ineq_constraint!(optmodel, CSTR)

alg = NLoptAlg(:LD_MMA)
opts = NLoptOptions(maxeval = 500, maxtime = 60)

res = optimize(optmodel, alg, x0, options = opts)
println("optimized volume: ", res.minimum, " m³  (",
    round(res.minimum / o0 * 100; digits=1), "% of initial)")
println("constraints satisfied: ", all(CSTR(res.minimizer) .<= 1e-6))

model2 = updatemodel(params, res.minimizer)

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
