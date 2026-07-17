# S4.2 gb_mma benchmark, OLD (publication) stack — replicates
# paper-scripts/S4.2_minvolume_spaceframe/gb_mma.jl verbatim: same model,
# variables, constraint rows, Nonconvex+NLopt LD_MMA, maxeval=1000,
# maxtime=300, ftol_rel=1e-3. Dumps results + the variable-structure
# indices the new-stack replication needs (OUTDIR env var).
using Pkg
const PUBDIR = joinpath(homedir(),
    "Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology",
    "Publications/DiffAnalysis_Structures_2024/DiffAnalysis_2024_CODE")
Pkg.activate(PUBDIR; io=devnull)

using LinearAlgebra
using DelimitedFiles

const OUTDIR = get(ENV, "OUTDIR", mktempdir())
isdir(OUTDIR) || mkpath(OUTDIR)

# patched include (same machinery as bench_old.jl): drop plotting/IO lines
const PATCHDIR = mktempdir()
const DROP = [r"^\s*import AsapToolkit", r"\batk\.", r"^\s*display\(",
    r"^\s*save\(", r"set_theme!", r"^\s*geo = Geo\(", r"Observable\(",
    r"@lift", r"^\s*GHsave\("]
function include_patched(mod, path)
    dir = dirname(abspath(path))
    lines = split(read(path, String), '\n')
    patched = map(l -> any(p -> occursin(p, l), DROP) ? "# [patched] " * l : l, lines)
    src = replace(join(patched, '\n'), "@__DIR__" => repr(dir))
    f = joinpath(PATCHDIR, basename(path))
    write(f, src)
    cd(() -> Base.include(mod, f), dir)
end

mod = Module(:PubGBMMA)
Core.eval(mod, :(using DiffAnalysis, LinearAlgebra))
include_patched(mod, joinpath(PUBDIR, "paper-scripts", "S4.2_minvolume_spaceframe", "init_problem.jl"))

println("MODEL nodes=", mod.model.nNodes, " elements=", mod.model.nElements,
    " nvars=", length(mod.x_init), " nstressed=", length(mod.i_stressed_elements))

# dump the problem structure for the new-stack replication
writedlm(joinpath(OUTDIR, "itop.txt"), mod.gen.itop)
writedlm(joinpath(OUTDIR, "i_stressed.txt"), mod.i_stressed_elements)
writedlm(joinpath(OUTDIR, "is_support.txt"),
    Int[mod.model.nodes[i].id == :support for i in 1:mod.model.nNodes])
writedlm(joinpath(OUTDIR, "x_init.txt"), mod.x_init)
writedlm(joinpath(OUTDIR, "lb.txt"), mod.params.lb)
writedlm(joinpath(OUTDIR, "ub.txt"), mod.params.ub)

# the gb_mma.jl optimization, verbatim
results = Core.eval(mod, quote
    function obj(x, p)
        geo = GeometricProperties(x, p)
        return dot(geo.L, geo.A)
    end
    OBJ = x -> obj(x, params)

    function cstr(x, p, dmax, fmax)
        res = solve_truss(x, p)
        vertical_displacements = res.U[3:3:end]
        axial_stresses = axial_stress(res, p)
        return [
            (-vertical_displacements .- dmax);
            (axial_stresses[i_stressed_elements] .- fmax);
        ]
    end
    CSTR = x -> cstr(x, params, dmax, fy)

    o0, do0 = withgradient(OBJ, x_init)
    c0, dc0 = withjacobian(CSTR, x_init)
    @assert all(c0 .< 0)

    alg = NLoptAlg(:LD_MMA)
    ftol = parse(Float64, get(ENV, "FTOL_REL", "1e-3"))
    opts = NLoptOptions(maxeval = 1000, maxtime = 300, ftol_rel = ftol)
    res = constrained_optimization(params, OBJ, CSTR, alg, opts)
    (o0, res)
end)
o0, res = results

println("RESULT initial_volume | ", o0)
println("RESULT obj_opt | ", res.obj_opt)
println("RESULT wall_time | ", res.time)
println("RESULT n_iter | ", length(res.obj_history))
println("RESULT stop | ", res.opt_stop_type)
c_opt = Core.eval(mod, :(CSTR($(res.x_opt))))
println("RESULT max_constraint | ", maximum(c_opt))

writedlm(joinpath(OUTDIR, "old_obj_history.txt"), res.obj_history)
writedlm(joinpath(OUTDIR, "old_x_opt.txt"), res.x_opt)
cmax_hist = [maximum(c) for c in res.cstr_history]
writedlm(joinpath(OUTDIR, "old_cmax_history.txt"), cmax_hist)
println("DONE")
