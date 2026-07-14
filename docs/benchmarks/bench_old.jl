# OLD-stack benchmark: publication environment (Asap 0.2.1 + DiffAnalysis's
# vendored differentiable layer — the AsapOptim-of-the-paper) on the SAME
# S4.2 spaceframe, same all-element area design vector, same objectives.

using Pkg
const PUBDIR = joinpath(homedir(),
    "Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology",
    "Publications/DiffAnalysis_Structures_2024/DiffAnalysis_2024_CODE")
Pkg.activate(PUBDIR; io=devnull)

using LinearAlgebra
using Asap
using DelimitedFiles

# build the S4.2 model VERBATIM via the publication's own init script
# (patched include, as used for fixture generation)
const PATCHDIR = mktempdir()
const DROP = [r"^\s*import AsapToolkit", r"\batk\.", r"^\s*display\(", r"^\s*save\("]
function include_patched(mod, path)
    dir = dirname(abspath(path))
    lines = split(read(path, String), '\n')
    patched = map(l -> any(p -> occursin(p, l), DROP) ? "# [patched] " * l : l, lines)
    src = replace(join(patched, '\n'), "@__DIR__" => repr(dir))
    f = joinpath(PATCHDIR, basename(path))
    write(f, src)
    cd(() -> Base.include(mod, f), dir)
end

mod = Module(:PubBench)
Core.eval(mod, :(using DiffAnalysis, LinearAlgebra))
include_patched(mod, joinpath(PUBDIR, "paper-scripts", "S4.2_minvolume_spaceframe", "init_problem.jl"))
Core.eval(mod, :(Asap.solve!(model)))

model = mod.model
println("MODEL nodes=$(model.nNodes) elements=$(model.nElements) freedofs=$(length(model.freeDOFs))")

function bench(f; n=10)
    f()
    times = Float64[]
    allocs = Float64[]
    for _ in 1:n
        GC.gc()
        a0 = Base.gc_bytes()
        t = @elapsed f()
        push!(times, t * 1e3)
        push!(allocs, (Base.gc_bytes() - a0) / 1e6)
    end
    return median(times), median(allocs)
end
median(v) = sort(v)[cld(length(v), 2)]

# 1. legacy linear analysis (rebuilds S from COO triplets every solve)
t, a = bench(() -> Asap.solve!(model; reprocess=true))
println("METRIC solve! | $t | $a")

# 2. legacy assembly alone
t, a = bench(() -> Asap.create_S!(model))
println("METRIC assemble | $t | $a")

# 3-5. legacy differentiable path with the SAME all-element-area variables
results = Core.eval(mod, quote
    A0 = [el.section.A for el in model.elements]
    vars = [AreaVariable(el, el.section.A, 1e-5, 1.0) for el in model.elements]
    params = TrussOptParams(model, TrussVariable[vars...])
    x0 = copy(params.values)

    obj_c(x) = begin
        res = solve_truss(x, params)
        dot(res.U, params.P)                       # legacy compliance form
    end
    obj_v(x) = begin
        geo = GeometricProperties(x, params)
        dot(geo.A, geo.L)
    end
    (params, x0, obj_c, obj_v)
end)
params, x0, obj_c, obj_v = results

t, a = bench(() -> Core.eval(mod, :(solve_truss($x0, $params))))
println("METRIC solve_structure | $t | $a")

g_c = Core.eval(mod, :(Zygote.gradient($obj_c, $x0)[1]))
g_v = Core.eval(mod, :(Zygote.gradient($obj_v, $x0)[1]))
t, a = bench(() -> Core.eval(mod, :(Zygote.gradient($obj_c, $x0))); n=5)
println("METRIC grad_compliance | $t | $a")
t, a = bench(() -> Core.eval(mod, :(Zygote.gradient($obj_v, $x0))); n=5)
println("METRIC grad_volume | $t | $a")

writedlm(joinpath(@__DIR__, "grad_compliance_old.txt"), g_c)
writedlm(joinpath(@__DIR__, "grad_volume_old.txt"), g_v)
println("DONE")
