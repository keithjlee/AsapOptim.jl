using AsapOptim, Asap, AsapToolkit
using Zygote
using LinearSolve, LinearAlgebra

# frame Optimization
begin
    section = Section(
        20e-3,
        2e8,
        8e7,
        7.95e-4,
        9.2e-5,
        3.11e-6
    )
end

begin
    Lx = 25.
    Ly = 15.
    x = 1.25
    y = 1.45
    z = 2.5
end

n = 24

# generate
begin
    support = :corner
    support_type = :pinned

    # loads
    load = [0., 0., -20]

    x_positions = range(0, Lx, n)
    y_positions = range(0, Ly, n)

    dx = Lx / (n-1)
    dy = Ly / (n-1)

    xyz = Vector{Vector{Float64}}()
    Xmatrix = zeros(Float64, n, n)
    Ymatrix = zeros(Float64, n, n)
    igrid = zeros(Int64, n, n)

    index = 1
    for iy = 1:n
        for ix = 1:n

            x = x_positions[iy]
            y = y_positions[ix]

            igrid[ix, iy] = index
            index += 1

            push!(xyz, [x, y, 0.])
            Xmatrix[ix, iy] = x
            Ymatrix[ix, iy] = y

        end
    end

    if support == :corner
        support_indices = [igrid[1, 1], igrid[n, 1], igrid[1, n], igrid[n, n]]
    elseif support == :x
        support_indices = igrid[[1, n], :]
    elseif support == :y
        support_indices = igrid[:, [1, n]]
    else
        support_indices = [igrid[[1, n], :][:]; igrid[2:n-1, [1, n]][:]]
    end

    #make nodes
    nodes = [Node(pos, :free, :free) for pos in xyz]

    #make support nodes
    for node in nodes[support_indices]
        fixnode!(node, support_type)
        node.id = :support
    end

    #make elements
    elements = Vector{Element}()

    #horizontal elements
    for i = 1:n
        for j = 1:n-1
            index = [igrid[i,j], igrid[i,j+1]]
            push!(elements, Element(nodes, index, section))
        end
    end

    #vertical elements
    for j = 1:n
        for i = 1:n-1
            index = [igrid[i,j], igrid[i+1,j]]
            push!(elements, Element(nodes, index, section)) 
        end
    end

    #loads
    loads = [NodeForce(node, load) for node in nodes[:free]]

    #assemble
    model = Model(nodes, elements, loads)
    Asap.solve!(model)
    # geo = Geo(model)
end;

# design variables
begin
    # n = 30
    @assert n % 2 == 0

    imid = Int(n / 2)

    iparent = igrid[2:imid, 2:imid]

    ichild1 = reverse(igrid[2:imid, imid+1:end-1], dims = 2)
    factors1 = [-1., 1.]

    ichild2 = reverse(igrid[imid+1:end-1, 2:imid], dims = 1)
    factors2 = [1., -1.]

    ichild3 = reverse(igrid[imid+1:end-1, imid+1:end-1])
    factors3 = [-1., -1.]
end

begin
    # make variables
    vars = FrameVariable[]
    coupled_vars = FrameVariable[]

    fac = .9
    x = dx * fac / 2
    y = dy * fac / 2
    z = 1.5


    for i in eachindex(iparent)

        i0 = iparent[i]
        i1 = ichild1[i]
        i2 = ichild2[i]
        i3 = ichild3[i]

        # x
        push!(vars, SpatialVariable(i0, 0., -x, x, :X))
        target = last(vars)

        push!(coupled_vars, CoupledVariable(i1, target, factors1[1]))
        push!(coupled_vars, CoupledVariable(i2, target, factors2[1]))
        push!(coupled_vars, CoupledVariable(i3, target, factors3[1]))

        # y
        push!(vars, SpatialVariable(i0, 0., -y, y, :Y))
        target = last(vars)

        push!(coupled_vars, CoupledVariable(i1, target, factors1[2]))
        push!(coupled_vars, CoupledVariable(i2, target, factors2[2]))
        push!(coupled_vars, CoupledVariable(i3, target, factors3[2]))

        # z
        # push!(vars, SpatialVariable(model.nodes[i0], 0., -z, z, :Z))
        # target = last(vars)

        # push!(vars, CoupledVariable(model.nodes[i1], target))
        # push!(vars, CoupledVariable(model.nodes[i2], target))
        # push!(vars, CoupledVariable(model.nodes[i3], target))
    end

    append!(vars, coupled_vars)
    # x, l, u  = AsapOptim.process_variables!(vars)
end


# iactive = findall(model.nodes, :free)

# begin
#     vars = FrameVariable[]
#     for node in model.nodes[iactive]
#         push!(vars, SpatialVariable(node.nodeID, 0., -x, x, :X))
#         push!(vars, SpatialVariable(node.nodeID, 0., -y, y, :Y))
#         push!(vars, SpatialVariable(node.nodeID, 0., -z, z, :Z))
#     end
# end

# vars = FrameVariable[
#     [SpatialVariable(node.nodeID, 0., -x, x, :X) for node in model.nodes[iactive]];
#     [SpatialVariable(node.nodeID, 0., -y, y, :Y) for node in model.nodes[iactive]];
#     [SpatialVariable(node.nodeID, 2., 0., z, :Z) for node in model.nodes[iactive]];
#     # [AreaVariable(element, 1e-5, .025) for element in model.elements];
#     # [SectionVariable(element, 1e-6, 1e-3, :Ix) for element in model.elements];
#     ]
# AsapOptim.process_variables!(vars)
# params = FrameOptParams2(model, vars);
params = FrameOptParams(model, vars)

#objective function
function objective_function(x::Vector{Float64}, p::FrameOptParams)

    res = solve_frame(x, p)

    dot(res.U, p.P)
end

OBJ = x -> objective_function(x, params)
@time g = gradient(OBJ, params.values)[1]

using Nonconvex, NonconvexNLopt

F = TraceFunction(OBJ)

omodel = Nonconvex.Model(F)
addvar!(
    omodel,
    params.lb,
    params.ub
)

alg = NLoptAlg(:LD_LBFGS)
opts = NLoptOptions(
    maxeval = 500,
    maxtime = 360,
    ftol_rel = 1e-8,
    ftol_abs = 1e-12,
    xtol_rel = 1e-8,
    xtol_abs = 1e-8
)


@time res = optimize(
    omodel,
    alg,
    params.values,
    options = opts
)

@show length(F.trace)

using kjlMakie; set_theme!(kjl_light_mono)
model2 = updatemodel(params, res.minimizer)
geo2 = Geo(model2)

begin
    dfac = Observable(0.)

    pts = @lift(Point3.(geo.nodes .+ $dfac .* geo.disp))
    els = @lift($pts[$geo.indices_flat])

    pts2 = @lift(Point3.(geo2.nodes .+ $dfac .* geo2.disp))
    els2 = @lift($pts2[geo2.indices_flat])
end

#visualize
begin
    fig = Figure(
        backgroundcolor = :white
    )

    ax = Axis3(
        fig[1,1],
        aspect = :data
    )

    asapstyle!(ax; ground = true)

    linesegments!(
        els,
        color = (:black, .25)
    )

    linesegments!(
        els2,
        # linewidth = geo2.Ix ./ maximum(geo2.Ix) .* 2
        color = geo2.Ix,
        # colormap = white2blue
    )

    # text!(
    #     pts,
    #     text = nvals
    # )

    sl = Slider(
        fig[2,1],
        startvalue = 0,
        range = range(0, 100, 250)
    )

    on(sl.value) do val
        dfac[] = val
    end

    on(dfac) do _
        autolimits!(ax)
    end

    fig
end