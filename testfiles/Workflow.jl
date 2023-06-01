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
        push!(vars, AreaVariable(el, tube.A / 2, 0., 35e3))
    end

    # # individual area optimization for web
    for element in model.elements[:web]
        push!(vars, AreaVariable(element, 500., 0., 20000.))
    end
end

@time problem = TrussOptParams(model, vars);

vals = problem.values

fac_compliance = model.compliance
fac_force = variation(getindex.(getproperty.(model.elements, :forces), 2))

function obj(values::Vector{Float64}, p::TrussOptParams)
    
    res = solvetruss(values, p)
    
    sum(res.L)
end

@time o1 = obj(vals, problem)
@time g1 = Zygote.gradient(var -> obj(var, problem), vals)[1]



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
        # maxiter = 1,
        )

end


res = OptimResults(problem, sol);
res.losstrace

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

amax = maximum((res.valtrace[end])[problem.indexer.iAg])
begin
    model1 = res.model

    i = Observable(1)
    lwfac = Observable(5)
    crfac = Observable(0.5)

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
end

stresses = AsapOptim.σaxial(getindex.(getproperty.(model1.elements, :forces), 2), getproperty.(getproperty.(model1.elements, :section), :A))

begin
    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax0); hidespines!(ax0)

    linesegments!(e0,
        color = :white,
        linestyle = :dash,
        linewidth = 0.5)

    linesegments!(e1,
        color = ff,
        colorrange = fr,
        colormap = pink2blue,
        linewidth = lw,
        )

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
    sleep(.01)
end