# initialize
begin
    using Asap, AsapToolkit, AsapOptim
    using kjlMakie; set_theme!(kjl_light)
    using Enzyme, Zygote
    using LinearAlgebra
end

# meta parameters
begin
    sec = rand(allHSSRound()) # choose a random round HSS section
    println(sec.name)
    asap_section = toASAPframe(sec, Steel_Nmm.E, Steel_Nmm.G) # generate an aSAP section
end

# geometry
begin
    nx = 13
    dx = 1000.
    dy = 1500.
end

#generate model
begin
    w2d = Warren2D(nx ,dx, dy, asap_section; load = [0., -40e3, 0.])
    model = w2d.model
    geo = Geo(model)
end

#visualize
# begin
#     dfac = Observable(0.)
#     cfac = Observable(.75)
#     lfac = Observable(5.)
#     mfac = Observable(15.)

#     pts = @lift(Point2.(geo.nodes_xy .+ $dfac .* geo.disp_xy))
#     els = @lift(vcat([$pts[id] for id in geo.indices]...))
#     cr = @lift($cfac .* (-1, 1) .* geo.max_abs_force)
#     lw = @lift(geo.areas ./ geo.max_area .* $lfac)
# end

# begin
#     fig = Figure(backgroundcolor = :white)
#     ax = Axis(
#         fig[1,1],
#         aspect = DataAspect()
#     )

#     hidespines!(ax)
#     tickstoggle!(ax)

#     ls_els = linesegments!(
#         els,
#         color = geo.forces,
#         colorrange = cr,
#         linewidth = lw,
#         colormap = pink2blue
#     )

#     sc_nds = scatter!(
#         pts,
#         color = :white,
#         strokecolor = :black,
#         markersize = mfac
#     )


#     sl = Slider(
#         fig[2,1],
#         startvalue = 0.,
#         range = 0:100
#     )

#     on(sl.value) do val
#         dfac[] = val
#     end

#     on(dfac) do _
#         autolimits!(ax)
#     end

#     fig
# end

#make parameters
begin
    vars = Vector{TrussVariable}()
    for node in model.nodes[:topchord]
        push!(vars, SpatialVariable(node, 0., -dy/2, dy, :Y))
    end

    for element in model.elements
        push!(vars, AreaVariable(element, element.section.A, 10., 20e3))
    end

    params = TrussOptParams(model, vars)
    x0 = params.values
end

#Allocation function
begin
    @time y_alloc = alloc(x0, params)
    @time ∇y_alloc = Zygote.gradient(x -> alloc(x, params), x0)[1]

    @show y_alloc
end;

#implicit non allocating function
begin
    prob0 = TrussOptProblem(model)

    @time y_nonalloc = nonalloc(x0, prob0, params)

    prob = TrussOptProblem(model)
    prob_collector = shadow(prob)

    # bx = zero(x0)

    # @time Enzyme.autodiff(
    #     Enzyme.Reverse, 
    #     nonalloc,
    #     Duplicated(x0, bx), 
    #     Duplicated(prob, prob_collector),
    #     Enzyme.Const(params)
    # )

    # ∇y_nonalloc = deepcopy(bx)

    @show y_nonalloc
end;

@show y_alloc - y_nonalloc

∇comp = [∇y_alloc ∇y_nonalloc]

prob_collector = shadow(prob)
dprob = Duplicated(prob, prob_collector)