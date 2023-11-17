using Asap, AsapToolkit, AsapOptim
using kjlMakie; set_theme!(kjl_light)
using Zygote, Nonconvex

Nonconvex.@load NLopt

# meta parameters
begin
    nx = 30 # number of bays in x
    dx = 1000. # x spacing
    ny = 10 # number of bays in 7
    dy = 1000. # y spacing
    dz = 1500. # spaceframe depth

    sec = rand(allHSSRound()) # choose a random round HSS section
    tube = toASAPtruss(sec, Steel_Nmm.E) # generate an aSAP section
end

# generate and analyze structural model
begin
    sf = SpaceFrame(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :x)
    model = sf.model
end

# define design variables
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
        push!(vars, AreaVariable(el, 5e3, 10., 20000.))
    end
end

# define params
params = TrussOptParams(model, vars)

# extract design variables
x0 = copy(params.values)

# define objective function
function obj_volume(x::Vector{Float64}, p::TrussOptParams)

    properties = GeometricProperties(x, p)

    dot(properties.A, properties.L)
end

o0 = obj_volume(x0, params)
@time g0 = Zygote.gradient(x -> obj_volume(x, params), x0)[1]

# define constraint
σ_max = 350.
d_max = min(nx * dx, ny * dy) / 240

function cstr_stress_disp(x::Vector{Float64}, p::TrussOptParams)

    #structural analysis
    res = solve_truss(x, p)

    #vertical displacements
    d_vertical = res.U[3:3:end]

    #axial stresses
    σ = axial_stress(res, p)

    #maximum vertical displacement
    d_critical = maximum(abs.(d_vertical))

    #maximum stress
    σ_critical = maximum(abs.(σ))

    #=
    Constraint g(x) in the form: g(x) + C ≤ 0
    =#

    return (
        d_critical - d_max,
        σ_critical - σ_max
    )

end

d_crit, sig_crit = cstr_stress_disp(x0, params)

# optimize
begin
    opt_alg = NLoptAlg(:LD_MMA)
    opt_options = NLoptOptions(
        maxeval = 250,
        ftol_rel = 1e-6,
        ftol_abs = 1e-6
    )
end

F = TraceFunction(x -> obj_volume(x, params))
opt_model = Nonconvex.Model(F)
addvar!(opt_model, params.lb, params.ub)

#solve
opt_results = Nonconvex.optimize(
    opt_model,
    opt_alg,
    x0,
    options = opt_options
)

losstrace = getproperty.(F.trace, :output)
xtrace = getproperty.(F.trace, :input)
xtrace_mat = hcat(xtrace...)
x_opt = opt_results.minimizer

o1 = obj_volume(x_opt, params)
c1 = cstr_stress_disp(x_opt, params)

