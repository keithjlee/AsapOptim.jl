# Enzyme FORWARD-mode probe — run in its own process: as of Enzyme 0.13.186
# a whole-pipeline forward pass can abort the process with an Enzyme compiler
# assertion (AdjointGenerator.h:318). The solve_free frule import itself is
# verified working (test_forward.jl); this probes the full pipeline.
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop([
    PackageSpec(path="/Users/keithlee/Documents/dev/Asap"),
    PackageSpec(path="/Users/keithlee/Documents/dev/AsapOptim"),
]; io=devnull)
Pkg.add([PackageSpec(name="Enzyme"), PackageSpec(name="ChainRulesCore"),
    PackageSpec(name="Zygote")]; io=devnull)
using Asap, AsapOptim, LinearAlgebra, ChainRulesCore
import Enzyme, Zygote

println("JULIA ", VERSION)

steel = Material(200e6, 80.0, 0.3)
section = Section(steel, 0.01)
n1 = Node([0.0, 0.0, 0.0], :pinned)
n2 = Node([4.0, 0.0, 0.0], vcat([true, true, false], trues(3)))
n3 = Node([8.0, 0.0, 0.0], :pinned)
n4 = Node([2.0, 3.0, 0.0], vcat([true, true, false], trues(3)), :top)
els = AbstractElement{Float64}[TrussElement(n1, n2, section), TrussElement(n2, n3, section),
    TrussElement(n1, n4, section), TrussElement(n4, n2, section)]
model = Asap.Model([n1, n2, n3, n4], els,
    AbstractLoad{Float64}[NodeForce(n2, [0.0, -50.0, 0.0])])
p = OptParams(model, AbstractVariable[SpatialVariable(n4, 0.0, -1.0, 1.5, :Y)])
x0 = copy(p.values)
obj = let p = p; x -> compliance(solve_structure(x, p), p); end

println("attempting whole-pipeline Enzyme forward JVP (may abort the process)...")
r = Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Forward), Enzyme.Const(obj),
    Enzyme.Duplicated, Enzyme.Duplicated(x0, [1.0]))
println("SUCCESS: JVP = ", r, "  zygote: ", Zygote.gradient(obj, x0)[1])
