using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra, FiniteDiff
using kjlMakie; set_theme!(kjl_dark)


### Create a spaceframe

#meta parameters
begin
    nx = 25
    dx = 750.
    ny = 25
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
# sf = generatespaceframe(nx, dx, ny, dy, dz, itp, tube, true; load = [0., 0., -30e3], support = :xy);

sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :y);
model = sf.truss;

#new random support locations
begin
    #remove existing supports
    for node in model.nodes[:support]
        fixnode!(node, :free)
        node.id = :bottom
    end

    #assign random supports
    for _ = 1:8
        iset = vec(rand(sf.isquares))

        for i in iset
            fixnode!(model.nodes[i], :pinned)
            model.nodes[i].id = :support
        end

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
    bottomAreaMaster = AreaVariable(model.elements[sf.ibottom[1]], 500., 0., 35000.)
    push!(vars, bottomAreaMaster)

    for i in sf.ibottom[2:end]
        push!(vars, CoupledVariable(model.elements[i], bottomAreaMaster))
    end

    # individual area optimization for web
    for element in model.elements[:web]
        push!(vars, AreaVariable(element, 500., 0., 20000.))
    end
end

@time problem = TrussOptParams(model, vars);

vals = problem.values
@time o1 = compliance(vals, problem);
@time g1 = Zygote.gradient(var -> compliance(var, problem), vals)[1];

func = Optimization.OptimizationFunction(compliance, Optimization.AutoZygote())
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

begin
    model1 = res1.model
    p1 = Point3.(getproperty.(model1.nodes, :position))
    e1 = vcat([p1[id] for id in getproperty.(model1.elements, :nodeIDs)]...)
    fo1 = getindex.(getproperty.(model1.elements, :forces), 2)
    a1 = getproperty.(getproperty.(model1.elements, :section), :A)

    a1normalized = a1 ./ maximum(a1)
end

begin
    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax0); hidespines!(ax0)

    linesegments!(e0,
        color = :white,
        linewidth = 2)

    ax1 = Axis3(fig[1,2],
        aspect = :data)

    hidedecorations!(ax1); hidespines!(ax1)

    linesegments!(e1,
        color = blue,
        linewidth = 5 .* a1normalized,
        )

    scatter!(p1[findall(model1.nodes, :support)],
        color = :white,
        markersize = 45)

    linkaxes!(ax0, ax1)

    axl = Axis(fig[2, 1:2],
        aspect = nothing,
        xlabel = "Iteration",
        ylabel = "Loss [Nmm]")

    lines!(res1.losstrace,
        color = :white,
        linewidth = 5)

    # ax2 = Axis(fig[2, 1:2],
    #     aspect = nothing)

    # series!(hcat(res1.valtrace...),
    #     solid_color = :white)

    fig
end