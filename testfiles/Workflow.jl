using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra
using kjlMakie; set_theme!(kjl_dark)


### Create a spaceframe

#meta parameters
begin
    nx = 15
    dx = 1000.
    ny = 15
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

# generate and extract model
sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :corner);
model = sf.truss;

# Optional random support generation
begin
    for node in model.nodes[:support]
        fixnode!(node, :free)
        node.id = :bottom
    end

    for _ = 1:10
        iset = vec(rand(sf.isquares))
        for node in model.nodes[iset]
            fixnode!(node, :pinned)
            node.id = :support
        end
    end

    for node in model.nodes[:support]
        node.position = node.position .- [0., 0., 5e3]
    end

    updateDOF!(model)
    solve!(model; reprocess = true)
end

#plot
begin
    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in getproperty.(model.elements, :nodeIDs)]...)

    f = getindex.(getproperty.(model.elements, :forces), 2)
    cr = maximum(abs.(f)) .* (-1, 1) .* .5

    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    linesegments!(e0, 
        color = f,
        colormap = pink2blue,
        colorrange = cr,
        linewidth = 5
        )

    scatter!(p0[findall(model.nodes, :support)],
        color = :yellow,
        strokecolor = :black,
        strokewidth = 5,
        markersize = 30)

    fig
end

#make variables
begin
    vars = Vector{TrussVariable}()

    # nodal Z variables
    lb = -50.
    ub = 3500.

    lbb = -1000.
    ubb = 500.

    #nodal xy variables for bottom nodes
    lbxy = -500.
    ubxy = 500.

    for node in model.nodes
        if node.id == :top
            push!(vars, SpatialVariable(node, 0., lb, ub, :Z))
            # push!(vars, SpatialVariable(node, 0., lbxy, ubxy, :X))
            # push!(vars, SpatialVariable(node, 0.,  lbxy, ubxy, :Y))
        end

        # if node.id == :bottom
        #     push!(vars, SpatialVariable(node, 0., lbb, ubb, :Z))
        # end

        # if node.id == :bottom
        #     push!(vars, SpatialVariable(node, 0., lbxy, ubxy, :X))
        #     push!(vars, SpatialVariable(node, 0.,  lbxy, ubxy, :Y))
        # end
    end

    #supports
    # for node in model.nodes[:support]
    #     push!(vars, SpatialVariable(node, 0., -6000., 250., :Z))
    # end

    # for el in model.elements[:bottom]
    #     push!(vars, AreaVariable(el, tube.A / 2, 0., 35e3))
    # end

    # # # individual area optimization for web
    # for element in model.elements[:web]
    #     push!(vars, AreaVariable(element, 500., 0., 20000.))
    # end

    for el in model.elements
        push!(vars, AreaVariable(el, 500., 0., 20000.))
    end

end

# generate a problem
problem = TrussOptParams(model, vars);

# extract design variables
vals = problem.values

# define objective function
function obj(values::Vector{Float64}, p::TrussOptParams)
    
    res = solvetruss(values, p)
    
    compliance(res, p)
end

# test objective
@time o1 = obj(vals, problem)
@time g1 = Zygote.gradient(var -> obj(var, problem), vals)[1];

#  define optimization problem
begin
    func = Optimization.OptimizationFunction(obj, 
        Optimization.AutoZygote(),
        )

    prob = Optimization.OptimizationProblem(func, vals, problem;
        lb = problem.lb,
        ub = problem.ub,
        )

    function cb(vals::Vector{Float64}, loss::Float64)
        push!(problem.losstrace, loss)
        push!(problem.valtrace, deepcopy(vals))
        false
    end

    cleartrace!(problem)
    @time sol = Optimization.solve(prob, 
        NLopt.LD_LBFGS();
        callback = cb,
        reltol = 1e-5,
        )

    # extract results
    res = OptimResults(problem, sol);
    res.losstrace
end


# intermediate values
begin
    fax = Vector{Vector{Float64}}()
    fmax = Vector{Float64}()

    for values in res.valtrace
        #populate values
        Xnew = AsapOptim.addvalues(problem.X, problem.indexer.iX, values[problem.indexer.iXg])
        Ynew = AsapOptim.addvalues(problem.Y, problem.indexer.iY, values[problem.indexer.iYg])
        Znew = AsapOptim.addvalues(problem.Z, problem.indexer.iZ, values[problem.indexer.iZg])
        Anew = AsapOptim.replacevalues(problem.A, problem.indexer.iA, values[problem.indexer.iAg])

        # vₑ
        elementvecs = AsapOptim.getevecs(Xnew, Ynew, Znew, problem)

        # Lₑ
        elementlengths = AsapOptim.getlengths(elementvecs)

        # vnₑ
        elementvecsnormalized = AsapOptim.getnormalizedevecs(elementvecs, elementlengths)

        # Γ
        rotmats = AsapOptim.getRmatrices(elementvecsnormalized)

        # kₑ
        klocs = AsapOptim.ktruss.(problem.E, Anew, elementlengths)

        # Kₑ
        kglobs = AsapOptim.getglobalks(rotmats, klocs)

        # K
        K = AsapOptim.assembleglobalK(kglobs, problem)

        # K⁻¹P
        disp = AsapOptim.solveU(K, problem)

        # U
        U = AsapOptim.replacevalues(zeros(problem.n), problem.freeids, disp)

        uE = [U[id] for id in problem.dofids]

        F = getindex.(rotmats .* kglobs .* uE, 2)

        push!(fax, F)
        push!(fmax, maximum(abs.(F)))
    end
end

# plotting assets
begin
    amax = maximum((res.valtrace[end])[problem.indexer.iAg])

    model1 = res.model

    i = Observable(1)
    lwfac = Observable(20)
    crfac = Observable(0.15)

    x = @lift(AsapOptim.addvalues(problem.X, problem.indexer.iX, (res.valtrace[$i])[problem.indexer.iXg]))
    y = @lift(AsapOptim.addvalues(problem.Y, problem.indexer.iY, (res.valtrace[$i])[problem.indexer.iYg]))
    z = @lift(AsapOptim.addvalues(problem.Z, problem.indexer.iZ, (res.valtrace[$i])[problem.indexer.iZg]))
    a = @lift(AsapOptim.replacevalues(problem.A, problem.indexer.iA, (res.valtrace[$i])[problem.indexer.iAg]))

    p1 = @lift(Point3.($x, $y, $z))
    e1 = @lift(vcat([$p1[id] for id in getproperty.(model1.elements, :nodeIDs)]...))

    lw = @lift($lwfac .* $a ./ amax)

    ff = @lift(fax[$i])
    fr = @lift($crfac .* (-1, 1) .* $fmax[$i])

    losstraceinc = @lift(res.losstrace[1:$i])

    supps = @lift($p1[findall(model.nodes, :support)])
end

begin
    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax0); hidespines!(ax0)

    init = linesegments!(e0,
        color = :white,
        linestyle = :dash,
        linewidth = 0.5)

    opt = linesegments!(e1,
        color = ff,
        colorrange = fr,
        colormap = pink2blue,
        linewidth = lw,
        )


    scatter!(supps,
        color = :yellow,
        strokecolor = :black,
        strokewidth = 5,
        markersize = 30)

    axl = Axis(fig[2, 1],
        aspect = nothing,
        xlabel = "Iteration",
        ylabel = "Loss [Nmm]")

    lines!(losstraceinc,
        color = :white,
        linewidth = 5)


    on(i) do _
        reset_limits!(ax0)
        reset_limits!(axl)
    end

    fig
end

for k in eachindex(res.losstrace)
    i[] = k
    sleep(.001)
end