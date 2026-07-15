#=
Example 4 — buying exactly as much connection stiffness as you need

USE CASE: a precast concrete frame. Beam-column moment connections are
made with grouted sleeves, corbels, or welded plates — and their rotational
stiffness is roughly proportional to fabrication cost (embedment length,
weld volume, bolt count). A fully rigid connection is expensive; a pinned
one is cheap but the frame drifts. The design question is therefore not
"rigid or pinned?" but: **what is the cheapest distribution of connection
stiffnesses that still meets the drift limit?**

That's a smooth optimization problem in the joint stiffnesses — exactly
what `JointVariable` + Asap's differentiable semi-rigid formulation
(Monforton–Wu end springs) makes tractable. Interior columns see different
demands than exterior ones, so the optimum is NOT uniform — the gradient
finds where stiffness actually earns its cost.

Self-contained: Asap, AsapOptim, Zygote only.
=#

using AsapOptim, Asap, LinearAlgebra
using Zygote

const VISUALIZE = get(ENV, "ASAP_EXAMPLE_VIZ", "true") == "true"
VISUALIZE && using CairoMakie

#=
A one-story, three-bay precast frame under lateral (wind) + gravity loads.
Beams connect to columns through semi-rigid joints at BOTH ends — and the
column bases are PINNED, so the beam-column connections are the frame's
ONLY lateral load path. Every kN·m/rad of drift resistance must be bought
at the joints; the drift limit genuinely binds.
=#
begin
    conc = Material(30e6, 2.4, 0.2)                    # kN, m — concrete
    col = Section(conc, 0.16, 2.13e-3, 2.13e-3, 3.6e-3)   # 400×400 column
    beam = Section(conc, 0.12, 1.6e-3, 0.9e-3, 2.0e-3)    # 300×400 beam

    H = 3.5          # story height
    B = 6.0          # bay width
    nbays = 3

    bases = [Node([B * (i - 1), 0.0, 0.0], :pinned, :base) for i in 1:nbays+1]
    tops = [Node([B * (i - 1), H, 0.0], :free, :top) for i in 1:nbays+1]

    columns = AbstractElement{Float64}[
        FrameElement(bases[i], tops[i], col, :column; rollangle = 0.0) for i in 1:nbays+1
    ]

    # beams: semi-rigid at both ends; k0 = a stiff starting guess
    k0 = 1e6                                            # kN·m/rad
    semirigid(k) = EndConditions(EndSprings(Inf, Inf, k, k), EndSprings(Inf, Inf, k, k))
    beams = AbstractElement{Float64}[
        FrameElement(tops[i], tops[i+1], beam, semirigid(k0), :beam; rollangle = 0.0) for i in 1:nbays
    ]

    # wind pushes the frame; gravity loads the beam-column joints
    # (nodal loads — element loads on jointed members would need their
    # fixed-end forces re-derived per design, which the pure path treats
    # as constant; AsapOptim guards this explicitly)
    W = 30.0                                            # lateral [kN]
    G = 120.0                                           # gravity per joint [kN]
    loads = AbstractLoad{Float64}[
        NodeForce(tops[1], [W, 0.0, 0.0]);
        [NodeForce(t, [0.0, -G, 0.0]) for t in tops]
    ]

    model = Asap.Model([bases; tops], [columns; beams], loads)
    planarize!(model)
    solve!(model)
end

drift_limit = H / 400
println("drift limit: ", round(drift_limit * 1000; digits=2), " mm")

#=
Design variables: one rotational stiffness per beam END (interior and
exterior joints may differ). The exterior pair is mirror-coupled; interior
joints share a group per side symmetry.
=#
begin
    klb, kub = 1e3, 1e8

    v_ext = JointVariable(beams[1], :start, k0, klb, kub)      # exterior joints
    v_int = JointVariable(beams[1], :end, k0, klb, kub)        # interior joints

    vars = AbstractVariable[
        v_ext,
        CoupledVariable((beams[end], :end), v_ext),            # mirror exterior
        v_int,
        CoupledVariable((beams[2], :both), v_int),             # center beam, both ends
        CoupledVariable((beams[end], :start), v_int),          # mirror interior
    ]

    params = OptParams(model, vars)
    x0 = copy(params.values)
end

# drift of the frame = max lateral displacement of the top nodes
top_ids = [n.index for n in model.nodes[:top]]
drift(res) = maximum(res.U[6*(i-1)+1] for i in top_ids)

#=
Objective: total connection stiffness (∝ fabrication cost) plus a smooth
exterior penalty on the drift limit. Two β continuation stages sharpen the
constraint. log-space design variables keep the search well-scaled across
five orders of magnitude of stiffness.
=#
function make_objective(params, β)
    return function (z)                    # z = log10(k)
        x = 10.0 .^ z
        res = solve_structure(x, params)
        cost = sum(x) / 1e6                # normalized stiffness budget
        viol = max(0.0, drift(res) / drift_limit - 1.0)
        return cost + β * viol^2
    end
end

function optimize_pgd(f, z0, lb, ub; iters = 80, α0 = 0.5, tol = 1e-10)
    z = clamp.(copy(z0), lb, ub)
    fz = f(z)
    for _ in 1:iters
        g = Zygote.gradient(f, z)[1]
        α = α0
        improved = false
        while α > 1e-12
            ztrial = clamp.(z .- α .* g, lb, ub)
            ftrial = f(ztrial)
            if ftrial < fz - tol * abs(fz)
                z, fz = ztrial, ftrial
                improved = true
                break
            end
            α /= 2
        end
        improved || break
    end
    return z, fz
end

zlb, zub = log10(klb) .* ones(2), log10(kub) .* ones(2)
z = log10.(x0)
for β in (1e2, 1e4)                        # constraint continuation
    global z, _ = optimize_pgd(make_objective(params, β), z, zlb, zub)
end
kopt = 10.0 .^ z

# evaluate the optimum
res_opt = solve_structure(kopt, params)
res_0 = solve_structure(x0, params)

println("\nJoint stiffness optimization (k in kN·m/rad)")
println("  uniform start:  exterior = interior = ", round(k0; sigdigits=3),
    "   drift = ", round(drift(res_0) * 1000; digits=2), " mm")
println("  optimized:      exterior = ", round(kopt[1]; sigdigits=3),
    ",  interior = ", round(kopt[2]; sigdigits=3))
println("                  drift = ", round(drift(res_opt) * 1000; digits=2),
    " mm  (limit ", round(drift_limit * 1000; digits=2), " mm)")
println("  stiffness budget: ", round(sum(kopt) / sum(x0) * 100; digits=1), "% of the uniform design")

@assert drift(res_opt) <= 1.02 * drift_limit "drift limit violated"

# materialize: the model now carries the optimized semi-rigid connections
model2 = updatemodel(params, kopt)
println("  beam 1 end conditions after update: ", model2.elements[nbays+2].ends)

if VISUALIZE
    # drift vs. stiffness-budget tradeoff curve, swept over uniform designs
    ks = 10.0 .^ range(3.5, 8, 30)
    drifts = [drift(solve_structure([k, k], params)) * 1000 for k in ks]

    fig = Figure()
    ax = Axis(fig[1, 1], xscale = log10,
        xlabel = "uniform joint stiffness [kN·m/rad]", ylabel = "drift [mm]")
    lines!(ks, drifts)
    hlines!([drift_limit * 1000], linestyle = :dash)
    scatter!([kopt[1], kopt[2]], fill(drift(res_opt) * 1000, 2), color = :red)
    display(fig)
end
