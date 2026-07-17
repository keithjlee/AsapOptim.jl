# S4.2 spaceframe (512 area vars), OLD registered stack (Asap 0.2.x + AsapOptim 0.1.3)
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.add([
    PackageSpec(name="AsapOptim", version="0.1"),
    PackageSpec(name="Zygote"),
    PackageSpec(name="BenchmarkTools"),
]; io=devnull)
Pkg.status(["Asap", "AsapOptim"])
using Asap, AsapOptim, LinearAlgebra, Zygote, BenchmarkTools

# rebuild the fixture on the OLD API (fixture dicts are package-agnostic)
include("/Users/keithlee/Documents/dev/Asap/test/characterization/fixtures_diffanalysis.jl")
def = DIFFANALYSIS_FIXTURES["S4.2_spaceframe/model"]
nodes = [TrussNode(Vector{Float64}(n["pos"]), Vector{Bool}(n["dof"])) for n in def["nodes"]]
els = [TrussElement(nodes[e["i"]], nodes[e["j"]],
    TrussSection(e["section"][1], e["section"][2])) for e in def["elements"]]
loads = [NodeForce(nodes[L["i"]], Vector{Float64}(L["value"])) for L in def["loads"]]
model = TrussModel(nodes, els, loads)
solve!(model)
println("MODEL nodes=", model.nNodes, " elements=", model.nElements,
    " freedofs=", length(model.freeDOFs))

# parity vs the pinned publication displacements
u_pin = Vector{Float64}(def["u"])
println("PARITY maxrel=", maximum(abs.(model.u .- u_pin)) / maximum(abs.(u_pin)))

vars = TrussVariable[AreaVariable(el, el.section.A, 1e-5, 1.0) for el in els]
params = TrussOptParams(model, vars)
x0 = copy(params.values)
OBJ = x -> begin res = solve_truss(x, params); dot(res.U, params.P) end
println("compliance = ", OBJ(x0))
b = @benchmark Zygote.withgradient($OBJ, $x0)
println("RESULT old-stack Zygote | ", round(median(b.times) / 1e6, digits=3),
    " ms | ", round(b.memory / 1e6, digits=2), " MB")
println("DONE")
