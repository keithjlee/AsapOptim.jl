using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra
using kjlMakie; set_theme!(kjl_dark)


### Create a spaceframe
#meta parameters
begin
    nx = 20
    dx = 1000.
    ny = 20
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

# generate and extract model
sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :xy);
model = sf.truss;

# generate random support locations and extend downwards

supportsquares = Vector{Vector{Int64}}()
begin
    for node in model.nodes[:support]
        fixnode!(node, :free)
        node.id = :bottom
    end

    for _ = 1:6
        iset = vec(rand(sf.isquares))
        push!(supportsquares, iset)
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

#make variables
begin
    vars = Vector{TrussVariable}()

    # nodal Z variables
    lb = -500.
    ub = 3500.

    lbb = -500.
    ubb = 500.

    #nodal xy variables for bottom nodes
    lbxy = -500.
    ubxy = 500.

    # spatial variables
    for node in model.nodes

        #top nodes can move in X, Y, Z
        if node.id == :top
            push!(vars, SpatialVariable(node, 0., lb, ub, :Z))
        end

        #bottom nodes can move in Z
        if node.id == :bottom
            push!(vars, SpatialVariable(node, 0., lbb, ubb, :Z))
        end
    end

    for el in model.elements
        push!(vars, AreaVariable(el, 500., 100., 20000.))
    end

    #each support square can move in XY
    lbxy2 = -3000.
    ubxy2 = 3000.

    for iset in supportsquares
        #master variable
        master_variable1 = SpatialVariable(model.nodes[iset[1]], 0., lbxy2, ubxy2, :X)
        master_variable2 = SpatialVariable(model.nodes[iset[1]], 0., lbxy2, ubxy2, :Y)

        
        push!(vars, master_variable1, master_variable2)
        n = length(vars)

        #tied variables
        for i in iset[2:end]
            push!(vars, CoupledVariable(model.nodes[i], vars[n-1]))
            push!(vars, CoupledVariable(model.nodes[i], vars[n]))
        end
    end

end

# generate a problem
problem = TrussOptParams(model, vars);

# extract design variables
vals = problem.values

# define objective function
function obj(values::Vector{Float64}, p::TrussOptParams)
    
    res = solvetruss(values, p)
    
    compliance(res, p) + dot(res.L, res.A) / 4
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
end


# extract results
res = OptimResults(problem, sol);
res.losstrace