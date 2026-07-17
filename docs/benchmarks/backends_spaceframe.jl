# S4.2 spaceframe (512 area vars), NEW v1.0 stack: Zygote / Mooncake / Enzyme
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop([
    PackageSpec(path="/Users/keithlee/Documents/dev/Asap"),
    PackageSpec(path="/Users/keithlee/Documents/dev/AsapOptim"),
]; io=devnull)
Pkg.add([PackageSpec(name="Zygote"), PackageSpec(name="Mooncake"),
    PackageSpec(name="Enzyme"), PackageSpec(name="DifferentiationInterface"),
    PackageSpec(name="BenchmarkTools"), PackageSpec(name="FiniteDifferences")]; io=devnull)

using Asap, AsapOptim, LinearAlgebra, BenchmarkTools
using DifferentiationInterface
import Zygote, Mooncake, Enzyme, FiniteDifferences

include("/Users/keithlee/Documents/dev/Asap/test/characterization/fixtures_diffanalysis.jl")
def = DIFFANALYSIS_FIXTURES["S4.2_spaceframe/model"]
nodes = [Node(Vector{Float64}(n["pos"]), vcat(Vector{Bool}(n["dof"]), trues(3))) for n in def["nodes"]]
els = [TrussElement(nodes[e["i"]], nodes[e["j"]],
    Section(Material(e["section"][2], 1.0, e["section"][3], 0.3), e["section"][1])) for e in def["elements"]]
loads = [NodeForce(nodes[L["i"]], Vector{Float64}(L["value"])) for L in def["loads"]]
model = Model(nodes, collect(AbstractElement{Float64}, els), collect(AbstractLoad{Float64}, loads))
solve!(model)
println("MODEL nodes=", length(nodes), " elements=", length(els),
    " freedofs=", length(model.cache.partition.free))

p = OptParams(model, AbstractVariable[AreaVariable(el, el.section.A, 1e-5, 1.0) for el in els])
x0 = copy(p.values)

make_obj(p) = x -> compliance(solve_structure(x, p), p)
obj = make_obj(p)

# reference: directional derivative vs FD (full 512-var FD grad is too slow)
v = normalize(cos.(1.0:length(x0)) .+ 0.1)
dref = FiniteDifferences.central_fdm(5, 1)(t -> obj(x0 .+ t .* v), 0.0)

backends = [
    ("Zygote", AutoZygote()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse),
        function_annotation=Enzyme.Const)),
]

for (name, be) in backends
    try
        prep = prepare_gradient(obj, be, x0)
        g = DifferentiationInterface.gradient(obj, prep, be, x0)
        err = abs(dot(g, v) - dref) / abs(dref)
        bb = @benchmark DifferentiationInterface.gradient($obj, $prep, $be, $x0)
        println("RESULT $name (DI, prepared) | ",
            round(median(bb.times) / 1e6, digits=3), " ms | ",
            round(bb.memory / 1e6, digits=2), " MB | dirderiv relerr=$err")
    catch e
        println("RESULT $name | FAIL: ", first(sprint(showerror, e), 300))
    end
end
println("DONE")
