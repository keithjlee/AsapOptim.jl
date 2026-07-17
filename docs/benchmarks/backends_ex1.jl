# example-1 truss, line-95 gradient: Zygote vs Mooncake vs Enzyme
# (DifferentiationInterface, prepared gradients), NEW v1.0 stack
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

# the example-1 Warren truss, exactly as in examples/truss-optimization1.jl
steel = Material(200e6, 80.0, 0.3)
section = Section(steel, 0.01)
L, nx, dy, load = 10.0, 12, 1.0, 75.0
dx = L / nx
bn = [Node([x, 0.0, 0.0], :free, :bottom) for x in range(0, L, nx + 1)]
fixnode!(first(bn), :pinned); fixnode!(last(bn), :xfree)
tn = [Node([dx / 2 + dx * (i - 1), dy, 0.0], :free, :top) for i in 1:nx]
els = AbstractElement{Float64}[
    [TrussElement(bn[i], bn[i+1], section, :bottom) for i in 1:nx];
    [TrussElement(tn[i], tn[i+1], section, :top) for i in 1:nx-1];
    [TrussElement(nb, nt, section, :web) for (nb, nt) in zip(bn[1:end-1], tn)];
    [TrussElement(nt, nb, section, :web) for (nt, nb) in zip(tn, bn[2:end])]]
loads = AbstractLoad{Float64}[NodeForce(n, [0.0, -load, 0.0]) for n in [bn; tn][:bottom]]
model = Asap.Model([bn; tn], els, loads)
planarize!(model); solve!(model)

xmin, xmax = (-1, 1) .* dx ./ 2 .* 0.9
ymin, ymax = -0.9dy, dy
vars = AbstractVariable[
    [SpatialVariable(n, 0.0, xmin, xmax, :X) for n in model.nodes[:top]];
    [SpatialVariable(n, 0.0, ymin, ymax, :Y) for n in model.nodes[:top]]]
p = TrussOptParams(model, vars)
x0 = copy(p.values)

make_obj(p) = x -> compliance(solve_structure(x, p), p)   # avoid Core.Box capture
obj = make_obj(p)
gref = FiniteDifferences.grad(FiniteDifferences.central_fdm(5, 1), obj, x0)[1]

# line 95 verbatim (unprepared Zygote), for continuity with the earlier table
b = @benchmark Zygote.withgradient($obj, $x0)
println("RESULT line95-verbatim Zygote.withgradient | ",
    round(median(b.times) / 1e3, digits=1), " μs | ",
    round(b.memory / 1e6, digits=2), " MB | maxrel=",
    maximum(abs.(Zygote.withgradient(obj, x0).grad[1] .- gref)) / maximum(abs.(gref)))

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
        err = maximum(abs.(g .- gref)) / maximum(abs.(gref))
        bb = @benchmark DifferentiationInterface.gradient($obj, $prep, $be, $x0)
        println("RESULT $name (DI, prepared) | ",
            round(median(bb.times) / 1e3, digits=1), " μs | ",
            round(bb.memory / 1e6, digits=2), " MB | maxrel=$err")
    catch e
        println("RESULT $name | FAIL: ", first(sprint(showerror, e), 300))
    end
end
println("DONE")
