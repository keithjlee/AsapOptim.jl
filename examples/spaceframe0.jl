using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra
using kjlMakie; set_theme!(kjl_dark)


### Create a spaceframe
#meta parameters
begin
    nx = 30
    dx = 1000.
    ny = 10
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

# generate and extract model
sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :x);
model = sf.truss;

#newloads
newloads = Vector{NodeForce}()
for node in model.nodes
    if node.position[2] >= ny*dy/2
        push!(newloads, NodeForce(node, [0., 0., -30e3]))
    end
end

Asap.solve!(model, [model.loads; newloads])

begin
    dfac = Observable(1.)
    p0 = @lift(Point3.(getproperty.(model.nodes, :position) .+ $dfac * getproperty.(model.nodes, :displacement)))
    e0 = @lift(vcat([$p0[id] for id in getproperty.(model.elements, :nodeIDs)]...))
    f0 = getindex.(getproperty.(model.elements, :forces), 2)
    c0 = maximum(abs.(f0)) .* (-1, 1) .* .2
    s0 = @lift($p0[findall(model.nodes, :support)])

    fig = Figure()
    ax = Axis3(fig[1,1],
        aspect = :data)

    gridtoggle!(ax); simplifyspines!(ax)

    linesegments!(e0,
        colormap = pink2blue,
        color = f0,
        colorrange = c0,
        linewidth = 4)

    scatter!(s0,
        color = :yellow,
        strokecolor = :black,
        strokewidth = 4,
        size = 100)
    

    fig
end

#make variables
begin
    vars = Vector{TrussVariable}()

    # nodal Z variables
    lb = -500.
    ub = 3500.

    lbb = -1000.
    ubb = 1500.

    #nodal xy variables for bottom nodes
    lbxy = -750.
    ubxy = 750.

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

    for el in model.elements
        push!(vars, AreaVariable(el, 250., 0., 20000.))
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
        reltol = 1e-3,
        )
end

# extract results
res = OptimResults(problem, sol);
res.losstrace

aset = [AsapOptim.replacevalues(problem.A, problem.indexer.iA, vals[problem.indexer.iAg]) for vals in res.valtrace]
anormalizer = maximum(maximum.(aset))

begin
    i = Observable(1)
    lfac = Observable(20)
    vals = @lift(res.valtrace[$i])

    x = @lift(AsapOptim.addvalues(problem.X, problem.indexer.iX, $vals[problem.indexer.iXg]))
    y = @lift(AsapOptim.addvalues(problem.Y, problem.indexer.iY, $vals[problem.indexer.iYg]))
    z = @lift(AsapOptim.addvalues(problem.Z, problem.indexer.iZ, $vals[problem.indexer.iZg]))
    a = @lift(AsapOptim.replacevalues(problem.A, problem.indexer.iA, $vals[problem.indexer.iAg]))

    lw = @lift($a ./ anormalizer .* $lfac)

    p = @lift(Point3.($x, $y, $z))
    e = @lift(vcat([$p[id] for id in problem.nodeids]...))
    s = @lift($p[findall(model.nodes, :support)])

end


begin
    fig = Figure(resolution = (1500,750))
    ax = Axis3(fig[1,1],
        aspect = :data)

    # gridtoggle!(ax); simplifyspines!(ax)
    hidedecorations!(ax); hidespines!(ax)

    linesegments!(e,
        color = :white,
        linewidth = lw)

    scatter!(s,
        color = :yellow,
        strokecolor = :black,
        strokewidth = 4,
        size = 100)
    

    on(i) do _
        reset_limits!(ax)
    end

    fig
end

for k = 1:res.niter
    i[] = k
    sleep(.05)
end

iterator = 1:res.niter

record(fig, "figures/canopy.gif", iterator; framerate = 20) do x
    i[] = x
end