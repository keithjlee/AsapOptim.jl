# initialize
begin
    using Asap, AsapToolkit, AsapOptim, LinearAlgebra, SparseArrays
    # using kjlMakie; set_theme!(kjl_light)
    using Enzyme, Zygote

    using LinearSolve
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
end;

#Allocation function

x0 = params.values ./ 5 .* (1 .+ rand())

begin
    @time y_alloc = alloc(x0, params)

    # @time dy_zygote = Zygote.gradient(x -> alloc(x, params), x0)[1]

    @show y_alloc
end;


#implicit non allocating function
begin

    prob_collector = shadow(prob)

    @time y_nonalloc = nonalloc(x0, prob, params)

    # bx1 = zero(x0)

    # @time Enzyme.autodiff(
    #     Enzyme.Reverse, 
    #     nonalloc,
    #     Duplicated(x0, bx1), 
    #     Duplicated(prob, prob_collector),
    #     Enzyme.Const(params)
    # )

    # dy_enzyme = deepcopy(bx1)

    @show y_nonalloc
end;

prob_collector = shadow(prob)

Ktest = deepcopy(prob.K[params.freeids, params.freeids])
@time obj = AsapOptim.solve_norm3(Ktest, params)

# THIS WORKS
dK = AsapOptim.explicit_zero(Ktest)
@time Enzyme.autodiff(
    Enzyme.Reverse,
    AsapOptim.solve_norm3,
    Duplicated(Ktest, dK),
    Enzyme.Const(params)
)


mutable struct TestStruct
    K::SparseMatrixCSC{Float64, Int64}
    Ksubset::SubArray{Float64, 2, SparseArrays.SparseMatrixCSC{Float64, Int64}, Tuple{Vector{Int64}, Vector{Int64}}, false}

    function TestStruct(model::Asap.AbstractModel)

        K = deepcopy(model.S)
        dK = @view(K[model.freeDOFs, model.freeDOFs])

        new(K, dK)
    end

    function TestStruct(Kref::SparseMatrixCSC{Float64, Int64}, indices::Vector{Int64})

        K = deepcopy(Kref)
        dK = @view(K[indices, indices])

        new(K, dK)
    end

end

import AsapOptim.shadow
function shadow(ts::TestStruct)
    K = deepcopy(ts.K)
    fill!(nonzeros(K), 0.)
    TestStruct(K, ts.Ksubset.indices[1])
end

function unorm(ts::TestStruct, params::TrussOptParams)

    lp = LinearProblem(ts.Ksubset, params.P[params.freeids])
    ls = LinearSolve.solve(lp, LUFactorization())

    norm(ls.u)
end

function F(x::Vector{Float64}, ts::TestStruct, params::TrussOptParams)

    ts.K.nzval .= x
    unorm(ts, params)

end

ts = TestStruct(model)
dts = shadow(ts)

x = deepcopy(ts.K.nzval)
dx = zero(x)

ts.K.nzval .= x

@time F(x, ts, params)

Enzyme.autodiff(
    Enzyme.Reverse,
    unorm,
    Duplicated(ts, dts),
    Enzyme.Const(params)
)