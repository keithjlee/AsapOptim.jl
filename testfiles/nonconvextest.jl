using Nonconvex

begin
    nx = 15
    dx = 750.
    ny = 15
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :xy);
model = sf.truss;

begin
    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in getproperty.(model.elements, :nodeIDs)]...)

    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    linesegments!(e0, color = :white)

    scatter!(p0[findall(model.nodes, :support)],
        markersize = 30)

    fig
end

#make variables
begin
    vars = Vector{TrussVariable}()

    # nodal Z variables
    lb = -250.
    ub = 2000.

    lbb = -4000.
    ubb = 250.

    #nodal xy variables for bottom nodes
    lbxy = -750.
    ubxy = 750.

    for node in model.nodes
        if node.id == :top
            push!(vars, SpatialVariable(node, 0., lb, ub, :Z))
        end

        if node.id == :bottom
            push!(vars, SpatialVariable(node, 0., lbb, ubb, :Z))
        end

        if node.id == :bottom
            push!(vars, SpatialVariable(node, 0., lbxy, ubxy, :X))
            push!(vars, SpatialVariable(node, 0.,  lbxy, ubxy, :Y))
        end
    end

    #supports
    for node in model.nodes[:support]
        push!(vars, SpatialVariable(node, 0., -6000., 250., :Z))
    end


    for el in model.elements[:bottom]
        push!(vars, AreaVariable(el, tube.A, 0., 35e3))
    end

    # # individual area optimization for web
    for element in model.elements[:web]
        push!(vars, AreaVariable(element, tube.A, 0., 20000.))
    end
end

@time problem = TrussOptParams(model, vars);

vals = problem.values

function objcompliance(values::Vector{Float64}, p::TrussOptParams)
    u = displacement(values, p)
    compliance(u, p)
end

#define model
optmodel = Nonconvex.Model(x -> objcompliance(x, problem))

#adding variables
addvar!(optmodel, problem.lb, problem.ub)

Nonconvex.@load NLopt

alg = NLoptAlg(:LD_LBFGS)
alg = NLoptAlg(:LD_MMA)
opts = NLoptOptions(ftol_rel = 1e-4)

@time r = Nonconvex.optimize(optmodel, alg, problem.values, options = opts)

model.compliance
r.minimum