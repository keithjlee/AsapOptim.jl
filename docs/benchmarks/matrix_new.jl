# AD comparison matrix, NEW stack (local Asap + AsapOptim):
# backend {Zygote, Mooncake, Enzyme-reverse, ForwardDiff}
#   × problem {example-1 truss (24 spatial vars), S4.2 spaceframe (512 area vars)}
#   × task {compliance gradient, full stress-constraint Jacobian}
# Prints RESULT lines: problem | task | backend | median ms | MB | maxrel vs Zygote.
# Run on each Julia version of interest. Enzyme FORWARD mode is probed by
# matrix_enzyme_forward_probe.jl (separate process — it can abort).
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop([
    PackageSpec(path="/Users/keithlee/Documents/dev/Asap"),
    PackageSpec(path="/Users/keithlee/Documents/dev/AsapOptim"),
]; io=devnull)
Pkg.add([PackageSpec(name="Zygote"), PackageSpec(name="Mooncake"),
    PackageSpec(name="Enzyme"), PackageSpec(name="ForwardDiff"),
    PackageSpec(name="DifferentiationInterface"),
    PackageSpec(name="BenchmarkTools")]; io=devnull)

using Asap, AsapOptim, LinearAlgebra, BenchmarkTools
using DifferentiationInterface
import Zygote, Mooncake, Enzyme, ForwardDiff

println("JULIA ", VERSION)

include("/Users/keithlee/Documents/dev/Asap/test/characterization/fixtures_diffanalysis.jl")

function ex1_problem()
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
    OptParams(model, vars)
end

function spaceframe_problem()
    def = DIFFANALYSIS_FIXTURES["S4.2_spaceframe/model"]
    nodes = [Node(Vector{Float64}(n["pos"]), vcat(Vector{Bool}(n["dof"]), trues(3))) for n in def["nodes"]]
    els = [TrussElement(nodes[e["i"]], nodes[e["j"]],
        Section(Material(e["section"][2], 1.0, e["section"][3], 0.3), e["section"][1])) for e in def["elements"]]
    loads = [NodeForce(nodes[L["i"]], Vector{Float64}(L["value"])) for L in def["loads"]]
    model = Model(nodes, collect(AbstractElement{Float64}, els), collect(AbstractLoad{Float64}, loads))
    solve!(model)
    OptParams(model, AbstractVariable[AreaVariable(el, el.section.A, 1e-5, 1.0) for el in els])
end

make_grad_obj(p) = x -> compliance(solve_structure(x, p), p)
make_jac_fn(p) = x -> axial_stress(solve_structure(x, p), p)

backends = [
    ("Zygote", AutoZygote()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("Enzyme-rev", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse),
        function_annotation=Enzyme.Const)),
    ("ForwardDiff", AutoForwardDiff()),
]

for (probname, p) in (("ex1-24svars", ex1_problem()), ("spaceframe-512avars", spaceframe_problem()))
    x0 = copy(p.values)
    obj = make_grad_obj(p)
    cst = make_jac_fn(p)
    g_ref = DifferentiationInterface.gradient(obj, AutoZygote(), x0)
    J_ref = DifferentiationInterface.jacobian(cst, AutoZygote(), x0)
    println("PROBLEM $probname: $(length(x0)) vars, $(length(J_ref[:, 1])) constraint rows")

    for (name, be) in backends
        # gradient of compliance
        try
            prep = prepare_gradient(obj, be, x0)
            g = DifferentiationInterface.gradient(obj, prep, be, x0)
            err = maximum(abs.(g .- g_ref)) / maximum(abs.(g_ref))
            b = @benchmark DifferentiationInterface.gradient($obj, $prep, $be, $x0)
            println("RESULT $probname | grad | $name | ",
                round(median(b.times) / 1e6, digits=4), " ms | ",
                round(b.memory / 1e6, digits=2), " MB | maxrel=$err")
        catch e
            println("RESULT $probname | grad | $name | FAIL: ", first(sprint(showerror, e), 150))
        end
        # full stress Jacobian
        try
            prep = prepare_jacobian(cst, be, x0)
            J = DifferentiationInterface.jacobian(cst, prep, be, x0)
            err = maximum(abs.(J .- J_ref)) / maximum(abs.(J_ref))
            b = @benchmark DifferentiationInterface.jacobian($cst, $prep, $be, $x0) seconds = 20
            println("RESULT $probname | jac  | $name | ",
                round(median(b.times) / 1e6, digits=4), " ms | ",
                round(b.memory / 1e6, digits=2), " MB | maxrel=$err")
        catch e
            println("RESULT $probname | jac  | $name | FAIL: ", first(sprint(showerror, e), 150))
        end
    end
end
println("DONE")
