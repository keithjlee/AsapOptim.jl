# NEW-stack benchmark: Asap v1.0 + AsapOptim v1.0 on the publication S4.2
# spaceframe (rebuilt exactly from the pinned fixture). Prints METRIC lines
# (name | median ms | allocated MB) and writes gradient vectors for
# cross-stack value comparison.

using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop([
        PackageSpec(path="/Users/keithlee/Documents/dev/Asap"),
        PackageSpec(path="/Users/keithlee/Documents/dev/AsapOptim"),
    ]; io=devnull)
Pkg.add(["Zygote", "FiniteDifferences"]; io=devnull)

using Asap, AsapOptim, LinearAlgebra, Zygote, FiniteDifferences
using DelimitedFiles

# ── rebuild the S4.2 spaceframe from the pinned publication fixture ─────────
include("/Users/keithlee/Documents/dev/Asap/test/characterization/fixtures_diffanalysis.jl")
def = DIFFANALYSIS_FIXTURES["S4.2_spaceframe/model"]

nodes = map(def["nodes"]) do n
    Node(Vector{Float64}(n["pos"]), vcat(Vector{Bool}(n["dof"]), trues(3)))
end
elements = map(def["elements"]) do e
    s = Vector{Float64}(e["section"])
    TrussElement(nodes[e["i"]], nodes[e["j"]], Section(Material(s[2], 1.0, s[3], 0.3), s[1]))
end
loads = map(def["loads"]) do L
    NodeForce(nodes[L["i"]], Vector{Float64}(L["value"]))
end
model = Model(nodes, collect(AbstractElement{Float64}, elements),
    collect(AbstractLoad{Float64}, loads))
solve!(model)
println("MODEL nodes=$(length(nodes)) elements=$(length(elements)) freedofs=$(length(model.cache.partition.free))")

# sanity: rebuilt model reproduces the pinned publication displacements
u_pin = Vector{Float64}(def["u"])
u_new = [model.results.u[6*(i-1)+c] for i in 1:length(nodes) for c in 1:3]
println("PARITY u_vs_publication maxrel=$(maximum(abs.(u_new .- u_pin) ./ (maximum(abs.(u_pin)))))")

# ── benchmark helpers ────────────────────────────────────────────────────────
function bench(f; n=10)
    f()                                    # warm
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

# 1. linear analysis (assembly + factorize + solve + post-process)
t, a = bench(() -> solve!(model))
println("METRIC solve! | $t | $a")

# 2. numeric re-assembly alone (the frozen-pattern path)
cache = model.cache
t, a = bench(() -> assemble_K!(cache, model))
println("METRIC assemble_K! | $t | $a")

# 3. differentiable forward solve (area design vector, all elements)
A0 = [el.section.A for el in elements]
vars = AbstractVariable[AreaVariable(el, el.section.A, 1e-5, 1.0) for el in elements]
p = OptParams(model, vars)
x0 = copy(p.values)
t, a = bench(() -> solve_structure(x0, p))
println("METRIC solve_structure | $t | $a")

# 4. gradients
obj_compliance(x) = compliance(solve_structure(x, p), p)
obj_volume(x) = begin
    geo = GeometricProperties(x, p)
    dot(geo.A, geo.L)
end

t, a = bench(() -> Zygote.gradient(obj_compliance, x0); n=5)
println("METRIC grad_compliance | $t | $a")
t, a = bench(() -> Zygote.gradient(obj_volume, x0); n=5)
println("METRIC grad_volume | $t | $a")

# 5. AD correctness: directional derivative vs central finite differences
g_c = Zygote.gradient(obj_compliance, x0)[1]
g_v = Zygote.gradient(obj_volume, x0)[1]
v = normalize(cos.(1.0:length(x0)) .+ 0.1)
fdm = central_fdm(5, 1)
dfd_c = fdm(t -> obj_compliance(x0 .+ t .* v), 0.0)
dfd_v = fdm(t -> obj_volume(x0 .+ t .* v), 0.0)
println("ADCHECK compliance zygote=$(dot(g_c, v)) fd=$dfd_c relerr=$(abs(dot(g_c,v)-dfd_c)/abs(dfd_c))")
println("ADCHECK volume zygote=$(dot(g_v, v)) fd=$dfd_v relerr=$(abs(dot(g_v,v)-dfd_v)/abs(dfd_v))")

writedlm(joinpath(@__DIR__, "grad_compliance_new.txt"), g_c)
writedlm(joinpath(@__DIR__, "grad_volume_new.txt"), g_v)
println("DONE")
