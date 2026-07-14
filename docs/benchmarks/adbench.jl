# Three-backend gradient benchmark: Zygote / Mooncake / Enzyme via
# DifferentiationInterface, on (a) a small mixed-variable truss MWE and
# (b) the publication S4.2 spaceframe (512 area variables).
using Pkg
Pkg.activate("adback"; io=devnull)
Pkg.resolve(io=devnull)
using Asap, AsapOptim, LinearAlgebra, SparseArrays
using DifferentiationInterface
import Zygote, Mooncake, Enzyme, FiniteDifferences

median(v) = sort(v)[cld(length(v), 2)]
function bench(f; n=7)
    f()
    ts = Float64[]; as = Float64[]
    for _ in 1:n
        GC.gc(); a0 = Base.gc_bytes(); t = @elapsed f()
        push!(ts, t * 1e3); push!(as, (Base.gc_bytes() - a0) / 1e6)
    end
    median(ts), median(as)
end

function small_problem()
    mat = Material(200e6, 1.0, 80.0, 0.3)
    sec = Section(mat, 1e-2)
    rot = [true, true, true]
    n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = Node([4.0, 0.0, 0.0], vcat([true, true, false], rot))
    n3 = Node([8.0, 0.0, 0.0], vcat([false, false, false], rot))
    n4 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot))
    n5 = Node([6.0, 3.0, 0.0], vcat([true, true, false], rot))
    els = AbstractElement{Float64}[
        TrussElement(n1, n2, sec), TrussElement(n2, n3, sec), TrussElement(n4, n5, sec),
        TrussElement(n1, n4, sec), TrussElement(n4, n2, sec), TrussElement(n2, n5, sec),
        TrussElement(n5, n3, sec)]
    loads = AbstractLoad{Float64}[NodeForce(n2, [0.0, -50.0, 0.0]), NodeForce(n5, [10.0, -30.0, 0.0])]
    model = Model([n1, n2, n3, n4, n5], els, loads)
    vars = AbstractVariable[
        SpatialVariable(n4, 0.0, -1.0, 1.5, :Y),
        AreaVariable(els[4], 1e-2, 1e-4, 5e-2),
        AreaVariable(els[1], 1e-2, 1e-4, 5e-2)]
    OptParams(model, vars)
end

function spaceframe_problem()
    include("/Users/keithlee/Documents/dev/Asap/test/characterization/fixtures_diffanalysis.jl")
    def = Main.DIFFANALYSIS_FIXTURES["S4.2_spaceframe/model"]
    nodes = [Node(Vector{Float64}(n["pos"]), vcat(Vector{Bool}(n["dof"]), trues(3))) for n in def["nodes"]]
    els = [TrussElement(nodes[e["i"]], nodes[e["j"]],
        Section(Material(e["section"][2], 1.0, e["section"][3], 0.3), e["section"][1])) for e in def["elements"]]
    loads = [NodeForce(nodes[L["i"]], Vector{Float64}(L["value"])) for L in def["loads"]]
    model = Model(nodes, collect(AbstractElement{Float64}, els), collect(AbstractLoad{Float64}, loads))
    solve!(model)
    OptParams(model, AbstractVariable[AreaVariable(el, el.section.A, 1e-5, 1.0) for el in els])
end

backends = [
    ("Zygote", AutoZygote()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse), function_annotation=Enzyme.Const)),
]
make_obj(p) = x -> compliance(solve_structure(x, p), p)   # avoids Core.Box capture
fdm = FiniteDifferences.central_fdm(5, 1)

for (probname, p) in (("small(3vars)", small_problem()), ("spaceframe(512vars)", spaceframe_problem()))
    x0 = copy(p.values)
    obj = make_obj(p)
    gref = FiniteDifferences.grad(fdm, obj, x0)[1]
    for (name, be) in backends
        try
            prep = prepare_gradient(obj, be, x0)
            g = DifferentiationInterface.gradient(obj, prep, be, x0)
            err = maximum(abs.(g .- gref)) / maximum(abs.(gref))
            t, a = bench(() -> DifferentiationInterface.gradient(obj, prep, be, x0))
            println("RESULT $probname | $name | $t ms | $a MB | maxrel=$err")
        catch e
            println("RESULT $probname | $name | FAIL: ", first(sprint(showerror, e), 200))
        end
    end
end
