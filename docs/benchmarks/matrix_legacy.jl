# AD comparison matrix, LEGACY stack (registered Asap 0.2.x + AsapOptim 0.1.3,
# Zygote only): compliance gradient + full stress Jacobian on both problems.
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.add([
    PackageSpec(name="AsapOptim", version="0.1"),
    PackageSpec(name="Zygote"),
    PackageSpec(name="BenchmarkTools"),
]; io=devnull)
using Asap, AsapOptim, LinearAlgebra, Zygote, BenchmarkTools

println("JULIA ", VERSION)
Pkg.status(["Asap", "AsapOptim"])

include("/Users/keithlee/Documents/dev/Asap/test/characterization/fixtures_diffanalysis.jl")

function ex1_problem()
    section = TrussSection(0.01, 200e6)
    L, nx, dy, load = 10.0, 12, 1.0, 75.0
    dx = L / nx
    bn = [TrussNode([x, 0.0, 0.0], :free, :bottom) for x in range(0, L, nx + 1)]
    fixnode!(first(bn), :pinned); fixnode!(last(bn), :xfree)
    tn = [TrussNode([dx / 2 + dx * (i - 1), dy, 0.0], :free, :top) for i in 1:nx]
    els = [
        [TrussElement(bn[i], bn[i+1], section, :bottom) for i in 1:nx];
        [TrussElement(tn[i], tn[i+1], section, :top) for i in 1:nx-1];
        [TrussElement(nb, nt, section, :web) for (nb, nt) in zip(bn[1:end-1], tn)];
        [TrussElement(nt, nb, section, :web) for (nt, nb) in zip(tn, bn[2:end])]]
    loads = [NodeForce(n, [0.0, -load, 0.0]) for n in [bn; tn][:bottom]]
    model = TrussModel([bn; tn], els, loads)
    planarize!(model); solve!(model)
    xmin, xmax = (-1, 1) .* dx ./ 2 .* 0.9
    ymin, ymax = -0.9dy, dy
    vars = TrussVariable[
        [SpatialVariable(n, 0.0, xmin, xmax, :X) for n in model.nodes[:top]];
        [SpatialVariable(n, 0.0, ymin, ymax, :Y) for n in model.nodes[:top]]]
    TrussOptParams(model, vars)
end

function spaceframe_problem()
    def = DIFFANALYSIS_FIXTURES["S4.2_spaceframe/model"]
    nodes = [TrussNode(Vector{Float64}(n["pos"]), Vector{Bool}(n["dof"])) for n in def["nodes"]]
    els = [TrussElement(nodes[e["i"]], nodes[e["j"]],
        TrussSection(e["section"][1], e["section"][2])) for e in def["elements"]]
    loads = [NodeForce(nodes[L["i"]], Vector{Float64}(L["value"])) for L in def["loads"]]
    model = TrussModel(nodes, els, loads)
    solve!(model)
    TrussOptParams(model, TrussVariable[AreaVariable(el, el.section.A, 1e-5, 1.0) for el in els])
end

for (probname, params) in (("ex1-24svars", ex1_problem()), ("spaceframe-512avars", spaceframe_problem()))
    x0 = copy(params.values)
    obj = x -> begin
        res = solve_truss(x, params)
        dot(res.U, params.P)
    end
    cst = x -> begin
        res = solve_truss(x, params)
        AsapOptim.axial_force(res, params) ./ res.A
    end
    obj(x0); cst(x0)   # warm

    b = @benchmark Zygote.withgradient($obj, $x0)
    println("RESULT $probname | grad | legacy-Zygote | ",
        round(median(b.times) / 1e6, digits=4), " ms | ",
        round(b.memory / 1e6, digits=2), " MB")

    Zygote.jacobian(cst, x0)
    b = @benchmark Zygote.jacobian($cst, $x0) seconds = 20
    println("RESULT $probname | jac  | legacy-Zygote | ",
        round(median(b.times) / 1e6, digits=4), " ms | ",
        round(b.memory / 1e6, digits=2), " MB")
end
println("DONE")
