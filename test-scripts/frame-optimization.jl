using AsapOptim, Asap, AsapToolkit
using Zygote
using LinearSolve, LinearAlgebra
using Profile

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

# generate
begin
    Lx = 25.
    Ly = 15.
    n = 28

    # loads
    load = [0., 0., -20]

    gridframe = GridFrame(Lx, n, Ly, n, section; load = load, support = :xy)
    model = gridframe.model
    geo = Geo(model)
end

# design variables

begin
    @assert n % 2 == 0

    imidx = Int(n / 2)
    imidy = Int(n / 2)

    iparent = gridframe.igrid[2:imidx, 2:imidy]

    ichild1 = reverse(gridframe.igrid[2:imidx, imidy+1:end-1], dims = 2)
    factors1 = [-1., 1.]

    ichild2 = reverse(gridframe.igrid[imidx+1:end-1, 2:imidy], dims = 1)
    factors2 = [1., -1.]

    ichild3 = reverse(gridframe.igrid[imidx+1:end-1, imidy+1:end-1])
    factors3 = [-1., -1.]
end

begin
    # make variables
    vars = FrameVariable[]

    fac = .9
    x = gridframe.dx * fac / 2
    y = gridframe.dy * fac / 2
    z = 1.5


    for i in eachindex(iparent)

        i0 = iparent[i]
        i1 = ichild1[i]
        i2 = ichild2[i]
        i3 = ichild3[i]

        # x
        push!(vars, SpatialVariable(model.nodes[i0], 0., -x, x, :X))
        target = last(vars)

        push!(vars, CoupledVariable(model.nodes[i1], target, factors1[1]))
        push!(vars, CoupledVariable(model.nodes[i2], target, factors2[1]))
        push!(vars, CoupledVariable(model.nodes[i3], target, factors3[1]))

        # y
        push!(vars, SpatialVariable(model.nodes[i0], 0., -y, y, :Y))
        target = last(vars)

        push!(vars, CoupledVariable(model.nodes[i1], target, factors1[2]))
        push!(vars, CoupledVariable(model.nodes[i2], target, factors2[2]))
        push!(vars, CoupledVariable(model.nodes[i3], target, factors3[2]))

        # z
        # push!(vars, SpatialVariable(model.nodes[i0], 0., -z, z, :Z))
        # target = last(vars)

        # push!(vars, CoupledVariable(model.nodes[i1], target))
        # push!(vars, CoupledVariable(model.nodes[i2], target))
        # push!(vars, CoupledVariable(model.nodes[i3], target))
    end
end

iactive = findall(model.nodes, :free)
vars = FrameVariable[
    # [SpatialVariable(node, 0., -1.25, 1.25, :X) for node in model.nodes[iactive]];
    # [SpatialVariable(node, 0., -1.25, 1.25, :Y) for node in model.nodes[iactive]];
    [SpatialVariable(node, 2., 0., 10., :Z) for node in model.nodes[iactive]];
    # [AreaVariable(element, 1e-5, .025) for element in model.elements];
    # [SectionVariable(element, 1e-6, 1e-3, :Ix) for element in model.elements];
    ]

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