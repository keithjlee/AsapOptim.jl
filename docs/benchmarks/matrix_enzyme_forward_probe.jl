# Enzyme FORWARD-mode probe — run in its own process AND under an external
# `timeout` (failure modes include ABORTS and compile-time HANGS, not
# exceptions). Status as of 2026-08-18, Enzyme 0.13.199 (see
# docs/AD_BACKENDS.md dated update):
#   * Julia 1.12.7: still ABORTS (AdjointGenerator.h:318 assertion),
#     identical to July.
#   * Julia 1.11.9 + default solver: works at ex1 scale (grad 0.30 ms,
#     Jacobian 0.31 ms, parity ~2e-14) but NEW REGRESSION at spaceframe
#     scale (512 vars): ReadOnlyMemoryError (0.13.186 ran it at 152 ms
#     in July).
#   * Re-pinning 0.13.186 on Julia 1.11.9 HANGS in Enzyme compilation of
#     the spaceframe forward Jacobian (>5 h CPU, killed) — the July
#     configuration no longer reproduces (environment drift).
# Status as of 2026-07-17 (Enzyme 0.13.186, native EnzymeRules.forward
# rule in AsapEnzymeExt), retained for history:
#   * Julia 1.11 + default solver: WORKS (machine-precision gradients and
#     Jacobians; ex1 Jacobian 0.30 ms, spaceframe 152 ms)
#   * Julia 1.11 + CachedSolver: ABORTS in Enzyme codegen around the
#     mutable solver struct (before rule dispatch) — use the default
#     solver with Enzyme, or ForwardDiff (which wins on speed anyway)
#   * Julia 1.12 (any solver): ABORTS (AdjointGenerator.h:318 assertion)
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop([
    PackageSpec(path="/Users/keithlee/Documents/dev/AsapEcosystem/Asap"),
    PackageSpec(path="/Users/keithlee/Documents/dev/AsapEcosystem/AsapOptim"),
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
