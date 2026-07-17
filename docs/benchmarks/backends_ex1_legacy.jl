# old stack, @benchmark median for consistency with backend table
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.add([
    PackageSpec(name="AsapOptim", version="0.1"),
    PackageSpec(name="Zygote"),
    PackageSpec(name="BenchmarkTools"),
]; io=devnull)
using Asap, AsapOptim, LinearAlgebra, Zygote, BenchmarkTools

A = 0.01; E = 200e6
section = TrussSection(A, E)
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
params = TrussOptParams(model, vars)
x0 = copy(params.values)
OBJ = x -> begin res = solve_truss(x, params); dot(res.U, params.P) end
b = @benchmark Zygote.withgradient($OBJ, $x0)
println("RESULT old-stack Zygote | ", round(median(b.times) / 1e3, digits=1),
    " μs median / ", round(minimum(b.times) / 1e3, digits=1), " μs min | ",
    round(b.memory / 1e6, digits=2), " MB")
