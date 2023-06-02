using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra
using kjlMakie; set_theme!(kjl_dark)


### Create a spaceframe
#meta parameters
begin
    nx = 20
    dx = 1000.
    ny = 15
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

# generate and extract model
sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :xy);
model = sf.truss;

#plot
begin
    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in getproperty.(model.elements, :nodeIDs)]...)

    f = getindex.(getproperty.(model.elements, :forces), 2)
    cr = maximum(abs.(f)) .* (-1, 1) .* .5

    fig0 = Figure(resolution = (2000,750))
    ax0 = Axis3(fig0[1,1],
        aspect = :data)

    ax0.ylabelvisible = ax0.xlabelvisible = ax0.zlabelvisible = false

    labelscale!(ax0, 1.5)

    simplifyspines!(ax0); gridtoggle!(ax0)

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

    fig0
end

#make variables
begin
    vars = Vector{TrussVariable}()

    # nodal Z variables
    lb = -500.
    ub = 3500.

    lbb = -1000.
    ubb = 500.

    #nodal xy variables for bottom nodes
    lbxy = -500.
    ubxy = 500.

    # spatial variables
    for node in model.nodes

        #top nodes can move in X, Y, Z
        if node.id == :top
            push!(vars, SpatialVariable(node, 0., lb, ub, :Z))
            push!(vars, SpatialVariable(node, 0., lbxy, ubxy, :X))
            push!(vars, SpatialVariable(node, 0.,  lbxy, ubxy, :Y))
        end

        #bottom nodes can move in Z
        if node.id == :bottom
            push!(vars, SpatialVariable(node, 0., lbb, ubb, :Z))
        end
    end

    #top element area variables
    for el in model.elements[:top]
        push!(vars, AreaVariable(el, 5000., 0., 20000.))
    end

    #bottom element area variables
    for el in model.elements[:bottom]
        push!(vars, AreaVariable(el, 5000., 0., 20000.))
    end

    #web elements must all have the same area
    iweb = findall(model.elements, :web)

    #define master variable
    web_variable_master = AreaVariable(model.elements[first(iweb)], 100., 0., 5000.)
    push!(vars, web_variable_master)

    #define coupled variables
    for i in iweb[2:end]
        push!(vars, CoupledVariable(model.elements[i], web_variable_master))
    end

end

# generate a problem
problem = TrussOptParams(model, vars);

# extract design variables
vals = problem.values

# define objective function
function obj(values::Vector{Float64}, p::TrussOptParams)
    
    res = solvetruss(values, p)
    
    compliance(res, p) # + 1e3 * sum(abs.(axialstress(res, p) .- 350))
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
        reltol = 1e-4,
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
    lwfac = Observable(5)
    crfac = Observable(0.25)

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

    trange = round.(range(0, res.solvetime, res.niter), digits = 3)
    ttext = @lift("Elapsed time: $(trange[$i])s")
end

begin
    fig = Figure(resolution = (2000, 1000))
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax0); hidespines!(ax0)

    # init = linesegments!(e0,
    #     color = :white,
    #     linestyle = :dash,
    #     linewidth = 0.5)

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
        markersize = 20)

    ax1 = Axis3(fig[1,2],
        elevation = 0,
        azimuth = pi/2,
        aspect = :data)

    hidedecorations!(ax1); hidespines!(ax1)

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
        markersize = 20)

    ax2 = Axis(fig[1,3],
        aspect = DataAspect())

    hidedecorations!(ax2); hidespines!(ax2)

    linesegments!(e1,
        color = ff,
        colorrange = fr,
        colormap = pink2blue,
        linewidth = lw,
        )

    scatter!(supps,
        color = :yellow,
        strokecolor = :black,
        strokewidth = 5,
        markersize = 20)

    axl = Axis(fig[2, 1:3],
        title = ttext,
        aspect = nothing,
        xlabel = "Iteration",
        ylabel = "Loss",
        yticklabelrotation = pi/4,
        ylabelpadding = 50,
        )

    labelscale!(axl, 1.5)

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
    sleep(.01)
end

iterator = collect(1:res.niter)

record(fig, "figures/spaceframe1_optim2.gif", iterator; framerate = 3) do x
    i[] = x
end