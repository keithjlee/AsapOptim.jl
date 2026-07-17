# S4.2 gb_mma benchmark, NEW v1.0 stack: identical problem (model from the
# pinned publication fixture, variable structure from gb_mma_old.jl's dumps
# in OUTDIR), identical MMA settings — but NLopt driven directly with
# reverse-mode (Zygote) objective, prepared FORWARD-mode (ForwardDiff)
# constraint Jacobian, and a shared CachedSolver factorization.
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop([
    PackageSpec(path="/Users/keithlee/Documents/dev/Asap"),
    PackageSpec(path="/Users/keithlee/Documents/dev/AsapOptim"),
]; io=devnull)
Pkg.add([PackageSpec(name="Zygote"), PackageSpec(name="ForwardDiff"),
    PackageSpec(name="DifferentiationInterface"), PackageSpec(name="NLopt")]; io=devnull)

using Asap, AsapOptim, LinearAlgebra, DelimitedFiles
using DifferentiationInterface
using NLopt
import Zygote, ForwardDiff

const OUTDIR = ENV["OUTDIR"]

# rebuild the exact publication model from the pinned fixture (node/element
# order preserved — displacement parity vs the publication is 1e-13)
include("/Users/keithlee/Documents/dev/Asap/test/characterization/fixtures_diffanalysis.jl")
def = DIFFANALYSIS_FIXTURES["S4.2_spaceframe/model"]
nodes = [Node(Vector{Float64}(n["pos"]), vcat(Vector{Bool}(n["dof"]), trues(3))) for n in def["nodes"]]
elements = [TrussElement(nodes[e["i"]], nodes[e["j"]],
    Section(Material(e["section"][2], 1.0, e["section"][3], 0.3), e["section"][1])) for e in def["elements"]]
loads = [NodeForce(nodes[L["i"]], Vector{Float64}(L["value"])) for L in def["loads"]]
model = Model(nodes, collect(AbstractElement{Float64}, elements),
    collect(AbstractLoad{Float64}, loads))
solve!(model)

# problem structure from the old run
itop = Int.(readdlm(joinpath(OUTDIR, "itop.txt")))
i_stressed = vec(Int.(readdlm(joinpath(OUTDIR, "i_stressed.txt"))))
is_support = vec(Bool.(readdlm(joinpath(OUTDIR, "is_support.txt"))))
x_init_old = vec(readdlm(joinpath(OUTDIR, "x_init.txt")))
lb_old = vec(readdlm(joinpath(OUTDIR, "lb.txt")))
ub_old = vec(readdlm(joinpath(OUTDIR, "ub.txt")))

# constants from init_problem.jl
n, d, dz = 8, 3.0, 2.25
fy = 350e3
dmax = n * d / 300
zmintop, zmaxtop = -0.9dz, dz
Amin, Amax = 1e-3, 0.2

# variables, replicating init_problem.jl exactly
vars = AbstractVariable[]
for i in diag(itop)
    is_support[i] && continue
    push!(vars, SpatialVariable(nodes[i], 0.0, zmintop, zmaxtop, :Z))
end
for i in 2:size(itop, 1), j in 1:i-1
    iparent, ichild = itop[i, j], itop[j, i]
    (is_support[iparent] || is_support[ichild]) && continue
    parent = SpatialVariable(nodes[iparent], 0.0, zmintop, zmaxtop, :Z)
    push!(vars, parent, CoupledVariable(nodes[ichild], parent))
end
for el in elements
    push!(vars, AreaVariable(el, 0.5Amax, Amin, Amax))
end

cs = Asap.CachedSolver()
params = OptParams(model, vars; solver = cs)
x0 = copy(params.values)

# parity with the publication problem: identical design vector and bounds
println("PARITY x_init  maxdiff = ", maximum(abs.(x0 .- x_init_old)))
println("PARITY bounds  maxdiff = ", max(maximum(abs.(params.lb .- lb_old)),
    maximum(abs.(params.ub .- ub_old))))

volume(x) = begin
    geo = GeometricProperties(x, params)
    dot(geo.L, geo.A)
end

# EXACT publication constraint rows (one-sided, unnormalized; v1.0 layout:
# vertical = Z = U[3:6:end])
constraints(x) = begin
    res = solve_truss(x, params)
    [
        (-res.U[3:6:end] .- dmax);
        (axial_stress(res, params)[i_stressed] .- fy)
    ]
end

o0 = volume(x0)
c0 = constraints(x0)
println("PARITY initial_volume = ", o0)
println("PARITY max_c0 = ", maximum(c0), "  feasible: ", all(c0 .< 0))

# NLopt callbacks: reverse objective, prepared forward Jacobian
obj_history = Float64[]
t_history = Float64[]
t_start = Ref(0.0)
function nlopt_objective(x::Vector, grad::Vector)
    if length(grad) > 0
        v, g = Zygote.withgradient(volume, x)
        grad .= g[1]
    else
        v = volume(x)
    end
    push!(obj_history, v)
    push!(t_history, time() - t_start[])
    return v
end

# JAC env var selects the constraint-Jacobian method:
#   "forwarddiff" (default) — prepared DifferentiationInterface ForwardDiff
#   "implicit"              — solution_tangents (implicit-function theorem)
const JAC = get(ENV, "JAC", "forwarddiff")
println("JACOBIAN method: ", JAC)

backend = AutoForwardDiff()
prep = prepare_jacobian(constraints, backend, x0)
function nlopt_constraints!(result::Vector, x::Vector, grad::Matrix)
    if length(grad) > 0
        if JAC == "implicit"
            t = solution_tangents(x, params)
            result .= [(-t.res.U[3:6:end] .- dmax);
                       (axial_stress(t.res, params)[i_stressed] .- fy)]
            grad .= vcat(-t.dU[3:6:end, :],
                         axial_stress_jacobian(t, params)[i_stressed, :])'
        else
            v, J = value_and_jacobian(constraints, prep, backend, x)
            result .= v
            grad .= J'
        end
    else
        result .= constraints(x)
    end
    return
end

m = length(c0)
opt = Opt(:LD_MMA, length(x0))
NLopt.lower_bounds!(opt, params.lb)
NLopt.upper_bounds!(opt, params.ub)
NLopt.min_objective!(opt, nlopt_objective)
NLopt.inequality_constraint!(opt, nlopt_constraints!, zeros(m))
NLopt.maxeval!(opt, 1000)
NLopt.maxtime!(opt, 300)
NLopt.ftol_rel!(opt, parse(Float64, get(ENV, "FTOL_REL", "1e-3")))

# warm up every code path OUTSIDE the timed run (the publication harness
# also differentiates once before optimizing)
Zygote.withgradient(volume, x0)
value_and_jacobian(constraints, prep, backend, x0)
if JAC == "implicit"
    t0 = solution_tangents(x0, params)
    axial_stress_jacobian(t0, params)
end

t_start[] = time()
wall = @elapsed begin
    minf, minx, ret = optimize(opt, x0)
    global minf, minx, ret
end

c_opt = constraints(minx)
println("RESULT initial_volume | ", o0)
println("RESULT obj_opt | ", minf)
println("RESULT wall_time | ", wall)
println("RESULT n_iter | ", NLopt.numevals(opt))
println("RESULT stop | ", ret)
println("RESULT max_constraint | ", maximum(c_opt))
println("RESULT cache | hits=", cs.hits, " refactorizations=", cs.refactorizations)

writedlm(joinpath(OUTDIR, "new_obj_history.txt"), obj_history)
writedlm(joinpath(OUTDIR, "new_t_history.txt"), t_history)
writedlm(joinpath(OUTDIR, "new_x_opt.txt"), minx)
println("DONE")
