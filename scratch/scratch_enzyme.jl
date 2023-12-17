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

    prob = TrussOptProblem(model; alg = KLUFactorization())

    OBJ_alloc = x -> alloc(x, params)
end;

#Allocation function

x0 = params.values ./ 5 .* (1 .+ rand())

begin
    @time y_alloc = alloc(x0, params)

    @time dy_zygote = Zygote.gradient(OBJ_alloc, x0)[1]

    @show y_alloc
end;


#implicit non allocating function
begin

    dprob = shadow(prob)

    @time y_nonalloc = nonalloc(x0, prob, params)

    bx1 = zero(x0)

    @time Enzyme.autodiff(
        Enzyme.Reverse, 
        nonalloc,
        Duplicated(x0, bx1), 
        Duplicated(prob, dprob),
        Enzyme.Const(params)
    )

    dy_enzyme = deepcopy(bx1)

    @show y_nonalloc
end;

@show norm(dy_zygote - dy_enzyme)

# step by step
function update_values!(x::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParams)

    #update values
    params.indexer.activeX && (prob.XYZ[params.indexer.iX, 1] += x[params.indexer.iXg])
    params.indexer.activeY && (prob.XYZ[params.indexer.iY, 2] += x[params.indexer.iYg])
    params.indexer.activeZ && (prob.XYZ[params.indexer.iZ, 3] += x[params.indexer.iZg])
    params.indexer.activeA && (prob.A[params.indexer.iA] = x[params.indexer.iAg])

    nothing

end

function update_element_vectors!(prob::TrussOptProblem, params::TrussOptParams)
    #element vectors
    prob.v = params.C * prob.XYZ
end

function update_element_properties!(prob::TrussOptProblem)

    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        prob.L[i] = norm(prob.v[i, :])
        prob.n[i, :] = prob.v[i, :] ./ prob.L[i]
    end

    nothing
end

function update_element_matrices!(prob::TrussOptProblem, params::TrussOptParams)
    #get local stiffness matrix, transformation matrix, and global stiffness matrix
    @simd for i in eachindex(prob.ke)
        prob.Γ[i][1, 1:3] = prob.n[i,:]
        prob.Γ[i][2, 4:6] = prob.n[i,:]
        prob.ke[i] = (params.E[i] * prob.A[i] / prob.L[i] * [1 -1; -1 1])
        prob.Ke[i] = prob.Γ[i]' * prob.ke[i] * prob.Γ[i]
    end

    nothing
end

prob_collector = shadow(prob)
x0 = params.values ./ 5 .* (1 .+ rand())

using JET
@report_opt update_values!(x0, prob, params)
@report_opt update_element_vectors!(prob, params)
@report_opt update_element_properties!(prob)
@report_opt update_element_matrices!(prob, params)
@report_opt AsapOptim.assemble_K!(prob, params)
@report_opt AsapOptim.solve_norm4(prob.K[params.freeids, params.freeids], params.P[params.freeids])

Ktest = deepcopy(prob.K[params.freeids, params.freeids])
@time obj = AsapOptim.solve_norm3(Ktest, params)

# THIS WORKS
@time dK = AsapOptim.explicit_zero(Ktest)
@time Enzyme.autodiff(
    Enzyme.Reverse,
    AsapOptim.solve_norm3,
    Duplicated(Ktest, dK),
    Enzyme.Const(params)
)

K = deepcopy(prob.K)
dK = AsapOptim.explicit_zero(K)

function sumsubset(K::SparseMatrixCSC{Float64, Int64}, inds::Vector{Int64})

    # norm(K[inds, inds])

    val = 0.
    for i in inds
        for j in inds
            val += K[i,j]
        end
    end

    val

end

@time sumsubset(K, params.freeids)

@time Enzyme.autodiff(Enzyme.Reverse, sumsubset, Enzyme.Active, Duplicated(K, dK), Enzyme.Const(params.freeids))

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