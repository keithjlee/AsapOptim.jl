using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra, FiniteDiff
using kjlMakie; set_theme!(kjl_dark)


### Create a spaceframe

#meta parameters
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

begin
    n = 5
    x = range(0, 1, n)
    y = range(0, 1, n)
    z = rand(n,n) .* 1500 .+ 500

    using Interpolations
    itp = cubic_spline_interpolation((x,y), z)
end
#generation and extraction
sf = generatespaceframe(nx, dx, ny, dy, dz, itp, tube, true; load = [0., 0., -30e3], support = :y);

sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :xy);
model = sf.truss;

begin
    for node in model.nodes[:support]
        fixnode!(node, :free)
        node.id = :bottom
    end

    for node in model.nodes[sf.iY1]
        fixnode!(node, :pinned)
        node.id = :support
    end

    for node in model.nodes[sf.iY2]
        fixnode!(node, :zfixed)
        node.id = :support
    end

    updateDOF!(model)
    solve!(model)
end

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
    lb = -1000.
    ub = 3500.

    lbb = -3500.
    ubb = 1000.

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

    # All bottom elements at once
    bottomAreaMaster = AreaVariable(first(model.elements[:bottom]), 500., 0., 35000.)
    push!(vars, bottomAreaMaster)

    for el in (model.elements[:bottom])[2:end]
        push!(vars, CoupledVariable(el, bottomAreaMaster))
    end

    # individual area optimization for web
    for element in model.elements[:web]
        push!(vars, AreaVariable(element, 500., 0., 20000.))
    end
end

@time problem = TrussOptParams(model, vars);

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

res = OptimResults(problem, sol);

amax = maximum((res.valtrace[end])[problem.indexer.iAg])
begin
    model1 = res.model

    i = Observable(1)
    lwfac = Observable(5)

    x = @lift(addvalues(problem.X, problem.indexer.iX, (res.valtrace[$i])[problem.indexer.iXg]))
    y = @lift(addvalues(problem.Y, problem.indexer.iY, (res.valtrace[$i])[problem.indexer.iYg]))
    z = @lift(addvalues(problem.Z, problem.indexer.iZ, (res.valtrace[$i])[problem.indexer.iZg]))
    a = @lift(replacevalues(problem.A, problem.indexer.iA, (res.valtrace[$i])[problem.indexer.iAg]))

    p1 = @lift(Point3.($x, $y, $z))
    e1 = @lift(vcat([$p1[id] for id in getproperty.(model1.elements, :nodeIDs)]...))

    lw = @lift($lwfac .* $a ./ amax)

    ztrace = stack([vals[problem.indexer.iZg] for vals in res.valtrace])
    ztraceinc = @lift(ztrace[:,1:$i])
    losstraceinc = @lift(res.losstrace[1:$i])

    atrace = stack([vals[problem.indexer.iAg] for vals in res.valtrace])
    atraceinc = @lift(atrace[:,1:$i])
end



begin
    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax0); hidespines!(ax0)

    linesegments!(e0,
        color = :white,
        linewidth = 2)


    # scatter!(p0,
    #     color = :white,
    #     markersize = 20)

    ax1 = Axis3(fig[1,2],
        aspect = :data)

    hidedecorations!(ax1); hidespines!(ax1)

    linesegments!(e1,
        color = blue,
        # linewidth = lw,
        )


    # scatter!(p1,
    #     color = blue,
    #     markersize = 20)

    linkaxes!(ax0, ax1)

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
        # reset_limits!(ax2)
    end

    fig
end

for k in eachindex(res.losstrace)
    i[] = k
    sleep(.01)
end