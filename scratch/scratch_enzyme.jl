# initialize
begin
    using Asap, AsapToolkit, AsapOptim
    using kjlMakie; set_theme!(kjl_light)
    using Enzyme, Zygote
    using LinearAlgebra, SparseArrays
end

# meta parameters
begin
    sec = rand(allHSSRound()) # choose a random round HSS section
    println(sec.name)
    asap_section = toASAPframe(sec, Steel_Nmm.E, Steel_Nmm.G) # generate an aSAP section
end

# geometry
begin
    nx = 15
    dx = 1000.

    ny = 25
    dy = 1500.

    dz = 2200.
end

#generate model
begin
    # gen = Warren2D(nx ,dx, dy, asap_section; load = [0., -40e3, 0.])
    gen = SpaceFrame(nx, dx, ny, dy, dz, asap_section; load = [0., 0., -10e3])

    model = gen.model
    geo = Geo(model)
end

#visualize
begin
    dfac = Observable(0.)
    cfac = Observable(.75)
    lfac = Observable(5.)
    mfac = Observable(15.)

    pts = @lift(Point2.(geo.nodes_xy .+ $dfac .* geo.disp_xy))
    els = @lift(vcat([$pts[id] for id in geo.indices]...))
    cr = @lift($cfac .* (-1, 1) .* geo.max_abs_force)
    lw = @lift(geo.areas ./ geo.max_area .* $lfac)
end

begin
    fig = Figure(backgroundcolor = :white)
    ax = Axis(
        fig[1,1],
        aspect = DataAspect()
    )

    hidespines!(ax)
    tickstoggle!(ax)

    ls_els = linesegments!(
        els,
        color = geo.forces,
        colorrange = cr,
        linewidth = lw,
        colormap = pink2blue
    )

    sc_nds = scatter!(
        pts,
        color = :white,
        strokecolor = :black,
        markersize = mfac
    )


    sl = Slider(
        fig[2,1],
        startvalue = 0.,
        range = 0:100
    )

    on(sl.value) do val
        dfac[] = val
    end

    on(dfac) do _
        autolimits!(ax)
    end

    fig
end

#make parameters
begin
    vars = Vector{TrussVariable}()
    # for node in model.nodes[:topchord]
    #     push!(vars, SpatialVariable(node, 0., -dy/2, dy, :Y))
    # end
    for node in model.nodes[:top]
        push!(vars, SpatialVariable(node, 0., -dz/2, dz, :Z))
    end

    for element in model.elements
        push!(vars, AreaVariable(element, element.section.A, 10., 20e3))
    end

    params = TrussOptParams(model, vars)
    x0 = params.values

    prob = TrussOptProblem(model)
end

#Allocation function
using BenchmarkTools
begin
    @time y_alloc = alloc(x0, params)

    @time dy_zygote = Zygote.gradient(x -> alloc(x, params), x0)[1]

    @show y_alloc
end;


#implicit non allocating function
begin

    prob = TrussOptProblem(model)
    prob_collector = shadow(prob)

    @time y_nonalloc = nonalloc(x0, prob, params)

    bx = zero(x0)

    @time Enzyme.autodiff(
        Enzyme.Reverse, 
        nonalloc,
        Duplicated(x0, bx), 
        Duplicated(prob, prob_collector),
        Enzyme.Const(params)
    )

    dy_enzyme = deepcopy(bx)

    @show y_nonalloc
end;

begin

    prob2 = TrussOptProblem2(model)

    @time y_nonalloc = nonalloc(x0, prob2, params)

    # bx = zero(x0)

    # @time Enzyme.autodiff(
    #     Enzyme.Reverse, 
    #     nonalloc,
    #     Duplicated(x0, bx), 
    #     Duplicated(prob, prob_collector),
    #     Enzyme.Const(params)
    # )

    # dy_enzyme = deepcopy(bx)

    @show y_nonalloc
end;

begin
    prob = TrussOptProblem(model)
    dprob = shadow(prob)
end

Enzyme.autodiff(
    Enzyme.Reverse,
    AsapOptim.assemble_K!,
    Enzyme.Const,
    Duplicated(prob, dprob),
    Enzyme.Const(params)
)