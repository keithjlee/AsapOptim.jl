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
    gen = SpaceFrame(nx, dx, ny, dy, dz, asap_section; load = [0., 0., -10e3])

    model = gen.model
    geo = Geo(model)
end

#make parameters
begin
    vars = Vector{TrussVariable}()
    
    for node in model.nodes[:top]
        push!(vars, SpatialVariable(node, 0., -dz/2, dz, :Z))
    end

    for element in model.elements
        push!(vars, AreaVariable(element, element.section.A, 10., 20e3))
    end

    params = TrussOptParamsNonalloc(model, vars)
    x0 = params.values

    params_alloc = TrussOptParams(model, vars)

    prob = TrussOptProblem(model)


    OBJ_alloc = x -> alloc(x, params_alloc)
end;

#Allocation function

x0 = params.values ./ 5 .* (1 .+ rand())

begin
    @time y_alloc = alloc(x0, params_alloc)

    @time dy_zygote = Zygote.gradient(OBJ_alloc, x0)[1]

    @show y_alloc
end;


#implicit non allocating function
begin

    prob = TrussOptProblem(model)
    @time y_nonalloc = nonalloc(x0, prob, params)

    dprob = shadow(prob)
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

#line by line
prob = TrussOptProblem(model)
dprob = shadow(prob)

#
AsapOptim.update_values!(x0, prob, params)

function f1(x0::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParamsNonalloc)

    AsapOptim.update_values!(x0, prob, params)

    norm(prob.A)
end

bx1 = zero(x0)
dprob = shadow(prob)
@time autodiff(
    Enzyme.Reverse,
    f1,
    Duplicated(x0, bx1),
    Duplicated(prob, dprob),
    Enzyme.Const(params)
)

@show bx1

function f2(x0::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParamsNonalloc)

    AsapOptim.update_values!(x0, prob, params)

    prob.v = params.C * prob.XYZ
    norm(prob.v)
end

bx2 = zero(x0)
dprob = shadow(prob)
@time autodiff(
    Enzyme.Reverse,
    f2,
    Duplicated(x0, bx2),
    Duplicated(prob, dprob),
    Enzyme.Const(params)
)

@show bx2

function f3(x0::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParamsNonalloc)

    AsapOptim.update_values!(x0, prob, params)

    prob.v = params.C * prob.XYZ
    
    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        @views prob.L[i] = norm(prob.v[i, :])
        @views prob.n[i, :] = prob.v[i, :] / prob.L[i]
    end

    dot(prob.L, prob.A)
end

bx3 = zero(x0)
dprob = shadow(prob)
@time autodiff(
    Enzyme.Reverse,
    f3,
    Duplicated(x0, bx3),
    Duplicated(prob, dprob),
    Enzyme.Const(params)
)

@show bx3

function f4(x0::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParamsNonalloc)

    AsapOptim.update_values!(x0, prob, params)

    prob.v = params.C * prob.XYZ
    
    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        @views prob.L[i] = norm(prob.v[i, :])
        @views prob.n[i, :] = prob.v[i, :] / prob.L[i]
    end

    #update element-wise stiffness/transformation information
    AsapOptim.element_update!(prob, params)

    sum(norm.(prob.Ke))
end

bx4 = zero(x0)
dprob = shadow(prob)
@time autodiff(
    Enzyme.Reverse,
    f4,
    Duplicated(x0, bx4),
    Duplicated(prob, dprob),
    Enzyme.Const(params)
)

@show bx4

function f5(x0::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParamsNonalloc)

    AsapOptim.update_values!(x0, prob, params)

    prob.v = params.C * prob.XYZ
    
    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        @views prob.L[i] = norm(prob.v[i, :])
        @views prob.n[i, :] = prob.v[i, :] / prob.L[i]
    end

    #update element-wise stiffness/transformation information
    AsapOptim.element_update!(prob, params)

    AsapOptim.assemble_K!(prob, params)

    norm(prob.K)
end

bx5 = zero(x0)
dprob = shadow(prob)
@time autodiff(
    Enzyme.Reverse,
    f5,
    Duplicated(x0, bx5),
    Duplicated(prob, dprob),
    Enzyme.Const(params)
)

@show bx5
@show dprob.Ke

function f6(x0::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParamsNonalloc)

    AsapOptim.update_values!(x0, prob, params)

    prob.v = params.C * prob.XYZ
    
    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        @views prob.L[i] = norm(prob.v[i, :])
        @views prob.n[i, :] = prob.v[i, :] / prob.L[i]
    end

    #update element-wise stiffness/transformation information
    AsapOptim.element_update!(prob, params)

    AsapOptim.assemble_K!(prob, params)

    LP = LinearProblem(copy(prob.K), prob.P)
    LS = LinearSolve.solve(LP)

    prob.u[params.freeids] = LS.u

    norm(prob.u)
end

bx6 = zero(x0)
prob = TrussOptProblem(model)
dprob = shadow(prob)

f6(x0, prob, params)

@time autodiff(
    Enzyme.Reverse,
    f6,
    Duplicated(x0, bx6),
    Duplicated(prob, dprob),
    Enzyme.Const(params)
)