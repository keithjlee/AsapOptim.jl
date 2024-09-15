using AsapOptim, LinearAlgebra

# You will need to also add the following packages for this example
using Asap, GLMakie, Nonconvex, NonconvexNLopt, Zygote

# initial section
begin
    A = 0.01 #m²
    E = 200E6 #KN/m²
    fy = 350e3 #kN/m²
    section = TrussSection(A, E)

    L = 10. # total span of truss
    nx = 12 # number of bays in span
    dy = 1.0 # truss depth
    load = 75. # force applied at each free node in bottom chord (downwards)
end

# Generate truss
begin
    dx = L / nx # bay length
    dmax = L / 360 #displacement limit

    # bottom nodes
    bottom_nodes = [TrussNode([x, 0., 0.], :free, :bottom) for x in range(0, L, nx+1)]

    # make supports
    fixnode!(first(bottom_nodes), :pinned)
    first(bottom_nodes).id = :pin

    fixnode!(last(bottom_nodes), :xfree)
    last(bottom_nodes).id = :roller

    # top nodes
    top_nodes = [TrussNode([dx/2 + dx * (i-1), dy, 0.], :free, :top) for i = 1:nx]

    # bottom chord elements
    bottom_elements = [TrussElement(bottom_nodes[i], bottom_nodes[i+1], section, :bottom) for i = 1:nx]

    #top chord elements
    top_elements = [TrussElement(top_nodes[i], top_nodes[i+1], section, :top) for i = 1:nx-1]

    # web elements
    web_elements = [
        [TrussElement(nb, nt, section, :web) for (nb, nt) in zip(bottom_nodes[1:end-1], top_nodes)];
        [TrussElement(nt, nb, section, :web) for (nt, nb) in zip(top_nodes, bottom_nodes[2:end])]
    ]

    # collect
    nodes = [bottom_nodes; top_nodes]
    elements = [bottom_elements; top_elements; web_elements]
    loads = [NodeForce(node, [0, -load, 0]) for node in nodes[:bottom]]

    # assemble
    model = TrussModel(nodes, elements, loads)
    planarize!(model)
    solve!(model)
end

# visualize
begin
    element_indices = vcat(Asap.nodeids.(model.elements)...)
    p = Point2.(getproperty.(model.nodes, :position))
    e = p[element_indices]

    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect())

    ylims!(-dy/2, 2dy)
    hidespines!(ax)
    hidedecorations!(ax)

    linesegments!(e, color = :black)
    scatter!(p, color = :white, strokecolor = :black, strokewidth = 1)

    fig
end

#=
CONSTRAINED MINIMUM VOLUME
=#

# make variables
xmin, xmax = (-1, 1) .* dx ./ 2 .* .9
ymin, ymax = -.9dy, dy

vars = TrussVariable[
    [SpatialVariable(node, 0., xmin, xmax, :X) for node in model.nodes[:top]];
    [SpatialVariable(node, 0., ymin, ymax, :Y) for node in model.nodes[:top]];
    [AreaVariable(element, element.section.A, .01element.section.A, 5element.section.A) for element in model.elements]
]

# make parameters
params = TrussOptParams(model, vars)
x0 = copy(params.values)

# objective function
function obj(x, p)
    geo = GeometricProperties(x, p)

    dot(geo.L, geo.A)
end

OBJ = x -> obj(x, params)
o0, ∇o0 = withgradient(OBJ, x0)

# constraint function
i_bottom_nodes = findall(model.nodes, :bottom)
function cstr(x, p, dmax, smax)
    res = solve_truss(x, p)
    
    vertical_displacements = res.U[2:3:end][i_bottom_nodes]
    stresses = AsapOptim.axial_force(res, p) ./ res.A

    return [
        (abs.(vertical_displacements) .- dmax);
        abs.(stresses) .- smax
    ]
end

CSTR = x -> cstr(x, params, dmax, fy)
c0, ∇c0 = withjacobian(CSTR, x0)

#=
Optimization
=#
F = TraceFunction(OBJ)
optmodel = Nonconvex.Model(F)
addvar!(optmodel, params.lb, params.ub)
add_ineq_constraint!(optmodel, CSTR)

alg = NLoptAlg(:LD_MMA)
opts = NLoptOptions(
    maxeval = 500,
    maxtime = 60
)

begin
    t0 = time()
    res = optimize(
        optmodel,
        alg,
        x0,
        options = opts
    )

    restime = time() - t0
end

# make new model from solution
model2 = updatemodel(params, res.minimizer)

# visualize
begin
    p2 = Point2.(getproperty.(model2.nodes, :position))
    e2 = p2[element_indices]

    areas = [e.section.A for e in model2.elements]
    lw = areas ./ maximum(areas) .* 4

    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect())

    ylims!(-dy/2, 3dy)
    hidespines!(ax)
    hidedecorations!(ax)

    # linesegments!(e, color = (:black, .1), linestyle = :dash)
    linesegments!(e2, color = :black, linewidth = lw)
    scatter!(p2, color = :white, strokecolor = :black, strokewidth = 1)

    fig
end

save("volume_solution.pdf", fig)

loss_history = getproperty.(F.trace, :output)

dt = range(0, restime, length(loss_history))

begin
    fig = Figure()
    ax = Axis(
        fig[1,1], 
        aspect = 3,
        xlabel = "TIME [s]",
        ylabel = "VOLUME [m³]",
        backgroundcolor = kjl_gray
    )

    graystyle!(ax)

    lines!(dt, loss_history)


    fig
end

save("volume_trace.pdf", fig)