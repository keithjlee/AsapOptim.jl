t2d = generatewarren2d(25, 2200., 3000., tube)

model = t2d.model

loads = [NodeForce(node, [0., -30e3, 0.]) for node in model.nodes[:bottomchord]]

solve!(model, loads)

begin

    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in getproperty.(model.elements, :nodeIDs)]...)

    fig = Figure()
    ax0 = Axis(fig[1,1],
        aspect = DataAspect())

    

    linesegments!(e0, color = :white)

    scatter!(p0[findall(model.nodes, :support)],
        markersize = 30)


    fig
end

#variables
begin
    vars = Vector{TrussVariable}()

    lbb = -4000.
    ubb = 1000.

    lbt = -1000.
    ubt = 2500.

    lbx = -1000.
    ubx = 1000.

    #bottom chord variabls
    for node in model.nodes[:bottomchord]
        push!(vars, SpatialVariable(node, 0., lbb, ubb, :Y))
        push!(vars, SpatialVariable(node, 0., lbx, ubx, :X))
    end

    #top chord variables
    for node in model.nodes[:topchord]
        push!(vars, SpatialVariable(node, 0., lbt, ubt, :Y))
    end

    #web areas
    awebMaster = AreaVariable(first(model.elements[:web]), tube.A, 500., 6000.)
    push!(vars, awebMaster)

    for el in (model.elements[:web])[2:end]
        push!(vars, CoupledVariable(el, awebMaster))
    end

end

problem = TrussOptParams(model, vars)

vals = problem.values

function obj(values::Vector{Float64}, p::TrussOptParams)
    u = displacement(values, p)
    compliance(u, p)
end

@time o1 = obj(vals, problem);
@time g1 = Zygote.gradient(var -> obj(var, problem), vals)[1];

func = Optimization.OptimizationFunction(obj, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(func, vals, problem;
    lb = problem.lb,
    ub = problem.ub)

function cb(vals::Vector{Float64}, loss::Float64)
    push!(problem.losstrace, loss)
    push!(problem.valtrace, deepcopy(vals))
    false
end

cleartrace!(problem)
@time sol = Optimization.solve(prob, NLopt.LD_LBFGS();
    callback = cb,
    reltol = 1e-4)

res1 = OptimResults(problem, sol);

amax = maximum((res1.valtrace[end])[problem.indexer.iAg])
begin
    model1 = res1.model

    i = Observable(1)
    lwfac = Observable(5)

    x = @lift(addvalues(problem.X, problem.indexer.iX, (res1.valtrace[$i])[problem.indexer.iXg]))
    y = @lift(addvalues(problem.Y, problem.indexer.iY, (res1.valtrace[$i])[problem.indexer.iYg]))
    z = @lift(addvalues(problem.Z, problem.indexer.iZ, (res1.valtrace[$i])[problem.indexer.iZg]))
    a = @lift(replacevalues(problem.A, problem.indexer.iA, (res1.valtrace[$i])[problem.indexer.iAg]))

    p1 = @lift(Point3.($x, $y, $z))
    e1 = @lift(vcat([$p1[id] for id in getproperty.(model1.elements, :nodeIDs)]...))

    lw = @lift($lwfac .* $a ./ amax)
end

ztrace = stack([vals[problem.indexer.iZg] for vals in res1.valtrace])
ztraceinc = @lift(ztrace[:,1:$i])
losstraceinc = @lift(res1.losstrace[1:$i])

atrace = stack([vals[problem.indexer.iAg] for vals in res1.valtrace])
atraceinc = @lift(atrace[:,1:$i])

begin
    fig = Figure()
    ax0 = Axis(fig[1,1],
        aspect = DataAspect())

    hidedecorations!(ax0); hidespines!(ax0)

    linesegments!(e0,
        color = :white,
        linewidth = 2)


    # scatter!(p0,
    #     color = :white,
    #     markersize = 20)

    ax1 = Axis(fig[1,2],
        aspect = DataAspect())

    hidedecorations!(ax1); hidespines!(ax1)

    linesegments!(e1,
        color = blue,
        linewidth = lw,
        )


    # scatter!(p1,
    #     color = blue,
    #     markersize = 20)

    # linkaxes!(ax0, ax1)

    axl = Axis(fig[2, 1:2],
        aspect = nothing,
        xlabel = "Iteration",
        ylabel = "Loss [Nmm]")

    lines!(losstraceinc,
        color = :white,
        linewidth = 5)

    # ax2 = Axis(fig[3, 1:2],
    #     aspect = nothing)

    # series!(atraceinc,
    #     solid_color = :white)


    on(i) do _
        reset_limits!(axl)
        reset_limits!(ax1)
    end

    fig
end