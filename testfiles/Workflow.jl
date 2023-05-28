using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra, FiniteDiff
using kjlMakie; set_theme!(kjl_dark)


### Create a spaceframe

#meta parameters
begin
    nx = 25
    dx = 750.
    ny = 15
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

n = 5
x = range(0, 1, n)
y = range(0, 1, n)
z = rand(n,n) .* 1500 .+ 500

using Interpolations
itp = cubic_spline_interpolation((x,y), z)

#generation and extraction
sf = generatespaceframe(nx, dx, ny, dy, dz, itp, tube, true; load = [0., 0., -30e3], support = :xy);

# sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :xy);
model = sf.truss;

newloads = Vector{NodeForce}()
for _ = 1:10
    iset = vec(rand(sf.isquares))

    for i in iset
        push!(newloads, NodeForce(model.nodes[i], [0., 0., -20e3]))
    end
end

solve!(model, [model.loads; newloads])

begin

    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in getproperty.(model.elements, :nodeIDs)]...)

    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    

    linesegments!(e0, color = :white)


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

function obj(values::Vector{Float64}, p::TrussOptParams)
    indexer = p.indexer

    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])

    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])

    ks = [kglobal(Xnew, Ynew, Znew, e, a, id) for (e, a, id) in zip(p.E, Anew, p.nodeids)]

    K = assembleglobalK(ks, p)

    U = solveU(K, p)

    return U' * p.P[p.freeids]
end

vals = problem.values
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

sol = Optimization.solve(prob, NLopt.LD_LBFGS();
    callback = cb,
    reltol = 1e-3)

res1 = OptimResults(problem, sol);

cleartrace!(problem)

sol2 = Optimization.solve(prob, NLopt.LD_MMA();
    callback = cb,
    reltol = 1e-3)

res2 = OptimResults(problem, sol2);

begin
    model1 = res1.model
    p1 = Point3.(getproperty.(model1.nodes, :position))
    e1 = vcat([p1[id] for id in getproperty.(model1.elements, :nodeIDs)]...)
    f1 = getindex.(getproperty.(model1.elements, :forces), 2)
    a1 = getproperty.(getproperty.(model1.elements, :section), :A)

    a1normalized = a1 ./ maximum(a1)

    model2 = res2.model
    p2 = Point3.(getproperty.(model2.nodes, :position))
    e2 = vcat([p2[id] for id in getproperty.(model2.elements, :nodeIDs)]...)
    f2 = getindex.(getproperty.(model2.elements, :forces), 2)
    a2 = getproperty.(getproperty.(model2.elements, :section), :A)

    a2normalized = a2 ./ maximum(a2)
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

    ax2 = Axis3(fig[1,3],
        aspect = :data)

    hidedecorations!(ax2); hidespines!(ax2)

    linesegments!(e2,
        color = green,
        linewidth = 5 .* a2normalized,
        )

    linkaxes!(ax0, [ax1, ax2])


    axloss = Axis(fig[2, 1:3],
        aspect = nothing)

    l1 = lines!(res1.losstrace,
        color = blue,
        linewidth = 4)

    l2 = lines!(res2.losstrace,
        color = green,
        linewidth = 4)

    fig
end



function obj2(values::Vector{Float64}, p::TrussOptProblem)
    indexer = p.indexer

    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])

    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])


    Ls = [L(Xnew[id]..., Ynew[id]..., Znew[id]...) for id in p.params.nodeids]
    kes = klocal.(problem.E, Anew, Ls)
    cxyz = [localvector(Xnew[id]..., Ynew[id]..., Znew[id]...) ./ l for (id, l) in zip(p.params.nodeids, Ls)]
    rs = [Rtruss(c...) for c in cxyz]

    ks = [r' * k * r for (r,k) in zip(rs, kes)]

    # ks = [kglobal(Xnew, Ynew, Znew, e, a, id) for (e, a, id) in zip(p.E, Anew, p.params.nodeids)]

    K = assembleglobalK(ks, p.params)

    U = solveU(K, p.params)

    U' * p.params.P[p.params.freeids]
end

@time o2 =  obj2(vals, problem);
@time g2 = Zygote.gradient(var -> obj2(var, problem), vals)[1];