#=
Example 3 — optimal geometry of a Vierendeel truss (FRAME optimization)

A Vierendeel truss has NO diagonals: it carries shear entirely through
frame action (bending of chords and posts through moment connections).
That makes it a genuinely frame-mechanical optimization problem — truss
idealizations can't even represent it (with pinned joints it would be a
mechanism).

We let the top-chord profile move vertically (mirror-symmetric) and ask:
what profile minimizes the volume-scaled compliance C·V — a size-invariant
stiffness metric — under uniform gravity loading?

Self-contained: needs only Asap, AsapOptim, Zygote. Optimization is a plain
projected gradient descent with backtracking, so you can see every moving
part. (Set ASAP_EXAMPLE_VIZ=true and add CairoMakie for the plots.)
=#

using AsapOptim, Asap, LinearAlgebra
using Zygote

const VISUALIZE = get(ENV, "ASAP_EXAMPLE_VIZ", "true") == "true"
VISUALIZE && using CairoMakie

#=
Build the Vierendeel: bottom chord at y = 0, top chord at y = h,
verticals ("posts") connecting them. Every member is a FrameElement with
rigid (:fixedfixed) joints — the load path IS the moment connections.
=#
begin
    steel = Material(200e6, 80.0, 0.3)                # kN, m
    chord = Section(steel, 6e-3, 8e-5, 8e-5, 1.6e-4)  # A, Ix, Iy, J
    post = Section(steel, 4e-3, 4e-5, 4e-5, 8e-5)

    L = 12.0        # span
    nbay = 6        # bays  (nbay + 1 posts)
    h = 1.5         # initial depth
    w = 40.0        # load per bottom node [kN]

    xs = collect(range(0, L, nbay + 1))
    bottom = [Node([x, 0.0, 0.0], :free, :bottom) for x in xs]
    top = [Node([x, h, 0.0], :free, :top) for x in xs]

    fixnode!(first(bottom), :pinned)
    fixnode!(last(bottom), :xfree)   # roller

    chords = AbstractElement{Float64}[
        [FrameElement(bottom[i], bottom[i+1], chord, :chord) for i in 1:nbay];
        [FrameElement(top[i], top[i+1], chord, :chord) for i in 1:nbay]
    ]
    posts = AbstractElement{Float64}[
        FrameElement(bottom[i], top[i], post, :post) for i in 1:nbay+1
    ]

    loads = AbstractLoad{Float64}[
        NodeForce(n, [0.0, -w, 0.0]) for n in bottom[2:end-1]
    ]

    model = Asap.Model([bottom; top], [chords; posts], loads)
    planarize!(model)
    solve!(model)
end

element_indices = vcat([[e.nodeStart.index, e.nodeEnd.index] for e in model.elements]...)

#=
Design variables: vertical position of each top node, ADDITIVE, coupled in
mirror-symmetric pairs so the optimizer preserves symmetry by construction.
=#
begin
    ymin, ymax = -0.75h, 2.0h
    npost = nbay + 1
    half = fld(npost, 2)

    vars = AbstractVariable[]
    for i in 1:half
        v = SpatialVariable(top[i], 0.0, ymin, ymax, :Y)
        push!(vars, v)
        push!(vars, CoupledVariable(top[npost+1-i], v))   # mirror partner
    end
    if isodd(npost)                                        # center post
        push!(vars, SpatialVariable(top[half+1], 0.0, ymin, ymax, :Y))
    end

    params = FrameOptParams(model, vars)   # alias of OptParams — the core is unified
    x0 = copy(params.values)
end

#=
Objective: C·V — compliance times volume. Pure compliance would just drive
every bound to its ceiling (more material is always stiffer); the product
rewards geometry that uses material WELL, the classic shape-optimization
trade.
=#
function objective(x, p)
    res = solve_structure(x, p)
    C = compliance(res, p)
    V = dot(res.A, res.L)
    return C * V
end
OBJ = x -> objective(x, params)

#=
Projected gradient descent with backtracking line search — deliberately
plain, so the whole optimization loop is visible.
=#
function optimize_pgd(f, x0, lb, ub; iters = 60, α0 = 1.0, tol = 1e-8)
    x = clamp.(copy(x0), lb, ub)
    fx = f(x)
    history = [fx]
    for it in 1:iters
        g = Zygote.gradient(f, x)[1]
        α = α0
        improved = false
        while α > 1e-10
            xtrial = clamp.(x .- α .* g, lb, ub)
            ftrial = f(xtrial)
            if ftrial < fx - tol * abs(fx)
                x, fx = xtrial, ftrial
                improved = true
                break
            end
            α /= 2
        end
        push!(history, fx)
        improved || break
    end
    return x, fx, history
end

@time f, dfdx = Zygote.withgradient(OBJ, x0)

f0 = OBJ(x0)
xopt, fopt, history = optimize_pgd(OBJ, x0, params.lb, params.ub)

println("Vierendeel geometry optimization")
println("  C·V initial:   ", round(f0; sigdigits=5))
println("  C·V optimized: ", round(fopt; sigdigits=5), "  (",
    round(fopt / f0 * 100; digits=1), "% of initial)")
println("  iterations:    ", length(history) - 1)

# materialize and inspect
model2 = updatemodel(params, xopt)
println("  top-chord profile [m]: ",
    round.([n.position[2] for n in model2.nodes[:top]]; digits=3))
println("  max deflection: ",
    round(minimum(displacement(model2.results, n)[2] for n in model2.nodes); sigdigits=4), " m")

if VISUALIZE
    p2 = Point2.([n.position[1:2] for n in model2.nodes])
    e2 = p2[element_indices]
    p1 = Point2.([n.position[1:2] for n in model.nodes])
    e1 = p1[element_indices]

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Vierendeel: initial (gray) vs optimized")
    hidespines!(ax); hidedecorations!(ax)
    linesegments!(e1, color = (:black, 0.2))
    linesegments!(e2, color = :black, linewidth = 2)
    scatter!(p2, color = :white, strokecolor = :black, strokewidth = 1)
    display(fig)

    fig2 = Figure()
    ax2 = Axis(fig2[1, 1], xlabel = "iteration", ylabel = "C·V")
    lines!(0:length(history)-1, history)
    display(fig2)
end
