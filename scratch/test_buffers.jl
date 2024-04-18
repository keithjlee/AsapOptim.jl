using AsapOptim

using AsapToolkit, Asap, LinearAlgebra, Zygote

#generate example structure
begin
    nx = 10
    dx = 1.25
    ny = 12
    dy = 1.5
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

x0 = copy(params.values) .+ rand(range(-dz, dz, 100), length(params.values))

z1 = AsapOptim.add_values(params.Z, params.indexer.iZ, x0[params.indexer.iZg])

z2 = AsapOptim.add_values_buffer(params.Z, params.indexer.iZ, x0[params.indexer.iZg])

@show norm(z1 - z2)

function func1(x::Vector{Float64}, p::TrussOptParams)

    z = AsapOptim.add_values(p.Z, p.indexer.iZ, x[p.indexer.iZg])

    return norm(z)
end

function func2(x::Vector{Float64}, p::TrussOptParams)

    z = AsapOptim.add_values_buffer(p.Z, p.indexer.iZ, x[p.indexer.iZg])

    return norm(z)
end

O1 = x -> func1(x, params)
O2 = x -> func2(x, params)

@time O1(x0)
@time O2(x0)

@time g1 = gradient(O1, x0)[1];
@time g2 = gradient(O2, x0)[1];

@show norm(g1 - g2)