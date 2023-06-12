using Asap, AsapToolkit, AsapOptim
using kjlMakie; set_theme!(kjl_dark)
import Nonconvex
Nonconvex.@load NLopt
Nonconvex.@load MMA
Nonconvex.@load Ipopt

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
    σ = 350.
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

#newloads
# newloads = Vector{NodeForce}()
# for node in model.nodes
#     if node.position[1] ≥ nx*dx/2 && node.position[2] ≥ ny*dy/2
#         push!(newloads, NodeForce(node, [0., 0., -40e3]))
#     end
# end

# iset = vec(rand(sf.isquares))
# L2 = [NodeForce(model.nodes[i], [0., 0., -100e3]) for i in iset]

# Asap.solve!(model, L2)

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
        push!(vars, AreaVariable(el, 10_000., 1_000., 20000.))
    end

end

# generate a problem
params = TrussOptParams(model, vars);

# define objective function
function obj(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Y = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Z = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    A = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    v = AsapOptim.getevecs(X, Y, Z, p)

    # Lₑ
    l = AsapOptim.getlengths(v)

    dot(A, l)
end

function cstr(values::Vector{Float64}, p::TrussOptParams)

    res = solvetruss(values, p)

    maximum(axialstress(res, p)) .- 350
end

# special structure to store traces
OBJ = x -> obj(x, params)
CSTR = x -> cstr(x, params)

#test
@time OBJ(params.values)
@time Zygote.gradient(OBJ, params.values)[1]

@time CSTR(params.values)
@time Zygote.gradient(CSTR, params.values)[1]

optmodel = Nonconvex.Model(OBJ)

#add variables and constraints
Nonconvex.addvar!(optmodel, params.lb, params.ub, init = copy(params.values))
Nonconvex.add_ineq_constraint!(optmodel, CSTR)

#define algorithm
begin
    alg = NLoptAlg(:LD_MMA)
    opts = NLoptOptions(ftol_rel = 1e-4)
    res = Nonconvex.optimize(optmodel, alg, params.values, options = opts)
end

begin
    alg2 = MMA()
    opts2 = MMAOptions()
    res2 = Nonconvex.optimize(optmodel, alg2, options = opts2)
end

begin
    alg3 = IpoptAlg()
    opts3 = IpoptOptions(first_order = true, tol = 1e-4)
    res3 = Nonconvex.optimize(optmodel, alg3, options = opts3)
end

#solution
begin
    sol = res.minimizer

    x = AsapOptim.addvalues(params.X, params.indexer.iX, sol[params.indexer.iXg])
    y = AsapOptim.addvalues(params.Y, params.indexer.iY, sol[params.indexer.iYg])
    z = AsapOptim.addvalues(params.Z, params.indexer.iZ, sol[params.indexer.iZg])
    a = AsapOptim.replacevalues(params.A, params.indexer.iA, sol[params.indexer.iAg])

    p1 = Point3.(x, y, z)
    e1 = vcat([p1[id] for id in params.nodeids]...)

    lfac = Observable(5.)
    lw = @lift(a ./ maximum(a) .* $lfac)
end

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

    oldpos = linesegments!(e0,
        colormap = pink2blue,
        color = f0,
        colorrange = c0,
        linewidth = 4)

    scatter!(s0,
        color = :yellow,
        strokecolor = :black,
        strokewidth = 4,
        size = 100)

    newpos = linesegments!(e1,
        color = :white,
        linewidth = lw)
    

    fig
end