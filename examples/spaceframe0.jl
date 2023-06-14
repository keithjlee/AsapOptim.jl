using Asap, AsapToolkit, AsapOptim

### Create a spaceframe
#meta parameters
begin
    nx = 25
    dx = 1000.
    ny = 15
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
end

#wavy
begin
    using Interpolations
    n = 5
    x = range(0, 1, n)
    y = range(0, 1, n)
    z = 3000 .* rand(n,n)

    itp = cubic_spline_interpolation((x,y), z)

    i = range(0,1, 50)
    j = range(0,1, 50)
    k = [itp(i,j) for i in i, j in j]
end

# generate and extract model
sf = generatespaceframe(nx, 
    dx, 
    ny, 
    dy, 
    dz,
    itp,
    tube,
    true; 
    load = [0., 0., -30e3], 
    support = :xy);

model = sf.truss;

#make variables
begin
    vars = Vector{TrussVariable}()

    # top node Z bounds
    lb = -500.
    ub = 3500.

    # bottom node Z bounds
    lbb = -1000.
    ubb = 1500.

    # bottom node XY bounds
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

    # area variables
    for el in model.elements
        push!(vars, AreaVariable(el, 20000., 0., 30000.))
    end

end

# generate a problem
params = TrussOptParams(model, vars);

# extract design variables
vals = params.values

# define objective function
function obj(values::Vector{Float64}, p::TrussOptParams)
    
    res = solvetruss(values, p)
    compliance(res, p)

end

# test objective
@time o1 = obj(vals, params)
@time g1 = Zygote.gradient(var -> obj(var, params), vals)[1]



#  define optimization problem
begin
    func = Optimization.OptimizationFunction(obj, 
        Optimization.AutoZygote(),
        )

    prob = Optimization.OptimizationProblem(func, vals, params;
        lb = params.lb,
        ub = params.ub,
        )

    function cb(vals::Vector{Float64}, loss::Float64)
        push!(params.losstrace, loss)
        push!(params.valtrace, deepcopy(vals))
        false
    end

    cleartrace!(params)
    @time sol = Optimization.solve(prob, 
        NLopt.LD_LBFGS();
        callback = cb,
        reltol = 1e-4,
        )
end

# extract results
res = OptimResults(params, sol)
@show res.losstrace