using AsapOptim

using AsapToolkit, Asap, LinearAlgebra, Zygote

#generate example structure
begin
    nx = 15
    dx = .75
    ny = 30
    dy = .5
    dz = .75

    A = 1e-3
    E = 200e6

    sec = TrussSection(A, E)

    load = [0., 0., -10.]

    sf = SpaceFrame(nx, dx, ny, dy, dz, sec; load = load)
    model = sf.model
end

# make variables
begin
    fac = .9
    vars = Vector{TrussVariable}()
    for node in model.nodes[:top]
        push!(vars, SpatialVariable(node, 0.,-fac * dz, dz, :Z))
    end
end

params = TrussOptParams(model, vars)

x0 = copy(params.values) #.+ rand(range(-dz, dz, 100), length(params.values))

function func(x::Vector{Float64}, p::TrussOptParams)

    res = AsapOptim.solve_truss_direct_buffer(x, p)

    return AsapOptim.compliance(res, p)
end

O = x -> func(x, params)

@time O(x0)

# @time g = gradient(O, x0)[1];

# using kjlMakie