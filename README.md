![](figures/gradients-axo.png)

# AsapOptim

This is a very much WIP package for high-performing *general* structural optimization in the [Asap.jl](https://github.com/keithjlee/Asap) environment. It provides the following:

- A set of data structures for controllable, extendable, and composable definitions of design variables and objective functions
- `CoupledVariable`s that reduce the dimensionality of design problems
- A complete set of custom chain rules for extremely high-performing automatic differentiation during gradient calculations
- A set of differentiable definitions of common structural objectives

# Small example
```julia
using Asap, AsapToolkit, AsapOptim
using Zygote, LinearAlgebra
using Nonconvex
Nonconvex.@load NLOpt


## Create a spaceframe

# meta parameters
nx = 30 # number of bays in x
dx = 1000. # x spacing
ny = 10 # number of bays in y
dy = 1000. # y spacing
dz = 1500. # spaceframe depth

sec = rand(allHSSRound()) # choose a random round HSS section
tube = toASAPtruss(sec, Steel_Nmm.E) # generate an aSAP section

# generate and extract model
sf = SpaceFrame(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :x)
model = sf.model

# Add a new set of loads on one half of structure
newloads = Vector{NodeForce}()
for node in model.nodes
    if node.position[2] >= ny*dy/2
        push!(newloads, NodeForce(node, [0., 0., -30e3]))
    end
end

# resolve with superimposed offset loads
Asap.solve!(model, [model.loads; newloads])

#make variables
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
for node in model.nodes[:top]
    push!(vars, SpatialVariable(node, 0., lb, ub, :Z))
    push!(vars, SpatialVariable(node, 0., lbxy, ubxy, :X))
    push!(vars, SpatialVariable(node, 0.,  lbxy, ubxy, :Y))
end

for node in model.nodes[:bottom]
    push!(vars, SpatialVariable(node, 0., lbb, ubb, :Z))
end

# area variables
for el in model.elements
    push!(vars, AreaVariable(el, 250., 1., 20000.))
end

# generate a problem
params = TrussOptParams(model, vars);

# extract design variables
x0 = params.values

# define objective function
function obj(values::Vector{Float64}, p::TrussOptParams)
    
    #get geometric properties
    properties = GeometricProperties(values, p)

    #total volume
    return dot(properties.A, properties.L)
end

# define constraints
σ_max = 350. #MPa
d_max = min(nx * dx, ny * dy) / 240

function cstr(values::Vector{Float64}, p::TrussOptParams)

    #structural analysis
    res = solve_truss(values, p)

    #vertical displacements
    d_vertical = res.U[3:3:end]

    #axial stresses
    σ = axial_stress(res, p)

    #maximum vertical displacement
    d_critical = minimum(d_vertical)

    #maximum stress
    σ_critical = maximum(abs.(σ))

    #=
    Constraint g(x) in the form: g(x) ≤ 0
    =#

    return (
        -d_critical - d_max,
        σ_critical - σ_max
    )
    
end

# define closures
OBJ = x -> obj(x, params)
CSTR = x -> cstr(x, params)

# define trace function
F = TraceFunction(OBJ)

# define an optimization problem
opt_model = Nonconvex.Model(F)

# add decision variables
addvar!(opt_model, params.lb, params.ub)

# add constraint
add_ineq_constraints!(opt_model, CSTR)

# define optimization options
opt_alg = NLoptAlg(:LD_MMA)

# define options
opt_options = NLoptOptions(
    maxeval = 100,
    ftol_rel = 1e-6,
    ftol_abs = 1e-6
)

#solve
opt_results = Nonconvex.optimize(
    opt_model,
    opt_alg,
    params.values,
    options = opt_options
)

```

