# AsapOptim v1.0 test suite: design-vector plumbing, parity of the
# optimization path with direct Asap solves, and end-to-end gradient checks
# through the core's pure path (Zygote vs central finite differences).

using AsapOptim
using Asap
using Test
using LinearAlgebra
using Zygote                # activates Asap's ChainRules extension
using ForwardDiff           # activates Asap's Dual solve path
using FiniteDifferences

const FDM = central_fdm(5, 1)

# a 2D-ish cantilevered truss with symmetric geometry variables + grouped areas
function testbed()
    mat = Material(200e6, 1.0, 80.0, 0.3)
    sec = Section(mat, 1e-2)
    rot = [true, true, true]

    n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = Node([4.0, 0.0, 0.0], vcat([true, true, false], rot))
    n3 = Node([8.0, 0.0, 0.0], vcat([false, false, false], rot))   # far pin: stable
    n4 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot), :top)
    n5 = Node([6.0, 3.0, 0.0], vcat([true, true, false], rot), :top)
    nodes = [n1, n2, n3, n4, n5]

    els = AbstractElement{Float64}[
        TrussElement(n1, n2, sec, :chord), TrussElement(n2, n3, sec, :chord),
        TrussElement(n4, n5, sec, :chord),
        TrussElement(n1, n4, sec, :web), TrussElement(n4, n2, sec, :web),
        TrussElement(n2, n5, sec, :web), TrussElement(n5, n3, sec, :web),
    ]
    loads = AbstractLoad{Float64}[
        NodeForce(n2, [0.0, -50.0, 0.0]), NodeForce(n5, [10.0, -30.0, 0.0])]
    model = Model(nodes, els, loads)

    # variables: top nodes move vertically (coupled symmetric), web areas grouped
    v_y = SpatialVariable(n4, 0.0, -1.0, 1.5, :Y)
    a_web = AreaVariable(els[4], 1e-2, 1e-4, 5e-2)
    vars = AbstractVariable[
        v_y,
        CoupledVariable(n5, v_y),                 # symmetric partner
        SpatialVariable(n4, 0.0, -0.5, 0.5, :X),
        a_web,
        CoupledVariable(els[5], a_web),           # grouped web areas
        CoupledVariable(els[6], a_web),
        CoupledVariable(els[7], a_web),
        AreaVariable(els[1], 1e-2, 1e-4, 5e-2),
    ]
    return model, vars
end

@testset "AsapOptim v1.0" begin
    model, vars = testbed()
    p = OptParams(model, vars)
    x0 = copy(p.values)

    @testset "parameter compilation" begin
        @test length(p.values) == 4                  # independents only
        @test p.lb <= p.values <= p.ub
        @test count(p.amask) == 5                    # 4 grouped webs + 1 chord
        @test TrussOptParams === OptParams === FrameOptParams
    end

    @testset "design evaluation ≡ direct Asap solve" begin
        # perturb the design, evaluate through the optimization path
        x = x0 .+ [0.4, -0.2, 2e-3, 5e-3]
        res = solve_structure(x, p)

        # materialize the same design and solve directly
        m2 = updatemodel(p, x)
        @test res.U ≈ m2.results.u rtol = 1e-10

        # node positions moved as declared (symmetric coupling honored)
        @test m2.nodes[4].position[2] ≈ 3.0 + 0.4
        @test m2.nodes[5].position[2] ≈ 3.0 + 0.4
        @test m2.nodes[4].position[1] ≈ 2.0 - 0.2
        @test m2.nodes[5].position[1] ≈ 6.0            # X var NOT coupled

        # grouped areas applied
        @test all(m2.elements[i].section.A ≈ x[3] for i in 4:7)
        @test m2.elements[1].section.A ≈ x[4]
        @test m2.elements[2].section.A ≈ 1e-2          # untouched

        # axial forces match Asap's recovered end forces
        Nopt = axial_force(res, p)
        for (i, el) in enumerate(m2.elements)
            @test Nopt[i] ≈ Asap.axial_force(m2.results, el) rtol = 1e-8 atol = 1e-8
        end

        # compliance consistent with the direct solve
        @test compliance(res, p) ≈ m2.results.compliance rtol = 1e-9
    end

    @testset "GeometricProperties (no solve)" begin
        x = x0 .+ [0.5, 0.0, 0.0, 0.0]
        geo = GeometricProperties(x, p)
        @test geo.L[4] ≈ norm([2.0, 3.5, 0.0])         # n1 → raised n4
        @test dot(geo.A, geo.L) > 0
    end

    @testset "gradients: volume objective (geometry only)" begin
        vol(x) = begin
            geo = GeometricProperties(x, p)
            dot(geo.A, geo.L)
        end
        g = Zygote.gradient(vol, x0)[1]
        gfd = FiniteDifferences.grad(FDM, vol, x0)[1]
        @test g ≈ gfd rtol = 1e-7
    end

    @testset "gradients: compliance through the solve" begin
        obj(x) = compliance(solve_structure(x, p), p)
        g = Zygote.gradient(obj, x0)[1]
        gfd = FiniteDifferences.grad(FDM, obj, x0)[1]
        @test g ≈ gfd rtol = 1e-6
        # more area ⇒ stiffer ⇒ lower compliance
        @test g[3] < 0 && g[4] < 0
    end

    @testset "gradients: stress-constraint style composite" begin
        obj(x) = begin
            res = solve_structure(x, p)
            σ = axial_stress(res, p)
            sum(abs2, σ) / 1e6
        end
        g = Zygote.gradient(obj, x0)[1]
        gfd = FiniteDifferences.grad(FDM, obj, x0)[1]
        @test g ≈ gfd rtol = 1e-6
    end

    @testset "frame model through the same path" begin
        mat = Material(200.0, 77.0, 1.0, 0.3)
        fsec = Section(mat, 1e4, 8e7, 3e7, 5e6)
        n1 = Node([0.0, 0.0, 0.0], :fixed)
        n2 = Node([3000.0, 0.0, 0.0], :free)
        n3 = Node([3000.0, 0.0, 3000.0], :fixed)
        beam = FrameElement(n1, n2, fsec, :beam)
        tie = TrussElement(n2, n3, Section(mat, 1e3), :tie)
        fmodel = Model([n1, n2, n3], AbstractElement{Float64}[beam, tie],
            AbstractLoad{Float64}[NodeForce(n2, [0.0, 0.0, -100.0])])

        fp = OptParams(fmodel, AbstractVariable[
            AreaVariable(tie, 1e3, 1e2, 1e4),
            SpatialVariable(n3, 0.0, -500.0, 500.0, :Z)])

        obj(x) = compliance(solve_frame(x, fp), fp)
        g = Zygote.gradient(obj, fp.values)[1]
        gfd = FiniteDifferences.grad(FDM, obj, fp.values)[1]
        @test g ≈ gfd rtol = 1e-6
        @test g[1] < 0                                  # thicker tie helps
    end
end

@testset "Network (FDM) optimization path" begin
    # a small hanging net: 3×3 grid, corners fixed
    ns = [Asap.FDMnode([Float64(i), Float64(j), 0.0], !(i in (0, 2) && j in (0, 2)))
          for j in 0:2 for i in 0:2]
    idx(i, j) = 3j + i + 1
    els = Asap.FDMelement[]
    for j in 0:2, i in 0:1
        push!(els, Asap.FDMelement(ns, idx(i, j), idx(i + 1, j), 1.0))
    end
    for j in 0:1, i in 0:2
        push!(els, Asap.FDMelement(ns, idx(i, j), idx(i, j + 1), 1.0))
    end
    loads = [Asap.FDMload(n, [0.0, 0.0, -1.0]) for n in ns if all(n.fixity)]
    network = Asap.Network(ns, els, loads)

    qv = QVariable(els[1], 1.5, 0.1, 10.0)
    vars = AbstractVariable[qv,
        CoupledVariable(els[2], qv),
        QVariable(els[7], 2.0, 0.1, 10.0)]
    np = NetworkOptParams(network, vars)
    @test length(np.values) == 2

    # forward parity with Asap's own FDM solver at the evaluated q
    x = [1.5, 2.0]
    res = solve_network(x, np)
    Asap.update_q!(network, collect(res.Q))
    Asap.solve!(network; reprocess = true)
    @test [res.X res.Y res.Z] ≈ vcat([n.position' for n in network.nodes]...) rtol = 1e-10
    @test all(isfinite, member_forces(res))

    # gradient of a smooth force-length objective
    obj(x) = begin
        r = solve_network(x, np)
        sum(abs2, r.Q .* r.L) / 100
    end
    g = Zygote.gradient(obj, x)[1]
    gfd = FiniteDifferences.grad(FDM, obj, x)[1]
    @test g ≈ gfd rtol = 1e-6
end

@testset "JointVariable (semi-rigid connection design)" begin
    mat = Material(200.0, 77.0, 1.0, 0.3)
    sec = Section(mat, 1e4, 8e7, 3e7, 5e6)
    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([3000.0, 0.0, 0.0], :free)
    n3 = Node([6000.0, 0.0, 0.0], :pinned)
    b1 = FrameElement(n1, n2, sec, EndConditions(EndSprings(Inf, Inf, 1e8, 1e8), rigid_end()), :b1; rollangle=0.0)
    b2 = FrameElement(n2, n3, sec, EndConditions(rigid_end(), EndSprings(Inf, Inf, 1e8, 1e8)), :b2; rollangle=0.0)
    model = Asap.Model([n1, n2, n3], AbstractElement{Float64}[b1, b2],
        AbstractLoad{Float64}[NodeForce(n2, [0.0, -100.0, 0.0])])

    jv = JointVariable(b1, :start, 1e8, 1e5, 1e12)
    vars = AbstractVariable[jv, CoupledVariable((b2, :end), jv)]   # mirrored pair
    p = OptParams(model, vars)
    x0 = copy(p.values)

    # evaluation matches a directly-built model at perturbed stiffness
    x = [3e8]
    res = solve_structure(x, p)
    m2 = updatemodel(p, x)
    @test res.U ≈ m2.results.u rtol = 1e-10
    @test m2.elements[1].ends.e1.kz ≈ 3e8
    @test m2.elements[2].ends.e2.kz ≈ 3e8    # coupled partner

    # physics: stiffer joints -> lower compliance
    c_soft = compliance(solve_structure([1e6], p), p)
    c_stiff = compliance(solve_structure([1e10], p), p)
    @test c_stiff < c_soft

    # gradient vs finite differences
    obj(x) = compliance(solve_structure(x, p), p)
    g = Zygote.gradient(obj, x0)[1]
    gfd = FiniteDifferences.grad(FDM, obj, x0)[1]
    # rtol reflects finite-difference truncation at the 1e8 stiffness scale
    @test g ≈ gfd rtol = 1e-4
    @test g[1] < 0                            # stiffer joint reduces compliance

    # the FEF-consistency guard fires for element loads on jointed members
    bn1 = Node([0.0, 0.0, 0.0], :fixed)
    bn2 = Node([3000.0, 0.0, 0.0], :pinned)
    bel = FrameElement(bn1, bn2, sec, :bad; rollangle=0.0)
    badmodel = Asap.Model([bn1, bn2], AbstractElement{Float64}[bel],
        AbstractLoad{Float64}[LineLoad(bel, [0.0, -1.0, 0.0])])
    @test_throws ErrorException OptParams(badmodel,
        AbstractVariable[JointVariable(bel, :start, 1e8, 1e5, 1e12)])
end

@testset "ForwardDiff forward mode" begin
    model, vars = testbed()
    p = OptParams(model, vars)
    x0 = copy(p.values)

    # compliance gradient: forward ≡ reverse to machine precision
    obj(x) = compliance(solve_structure(x, p), p)
    g_fwd = ForwardDiff.gradient(obj, x0)
    g_rev = Zygote.gradient(obj, x0)[1]
    @test g_fwd ≈ g_rev rtol = 1e-12

    # stress-constraint JACOBIAN — the forward-mode sweet spot
    # (n_outputs = n_elements ≥ n_variables)
    cstr(x) = axial_stress(solve_structure(x, p), p)
    J_fwd = ForwardDiff.jacobian(cstr, x0)
    J_rev = Zygote.jacobian(cstr, x0)[1]
    @test J_fwd ≈ J_rev rtol = 1e-12

    # volume objective through GeometricProperties (no solve)
    vol(x) = begin
        geo = GeometricProperties(x, p)
        dot(geo.A, geo.L)
    end
    @test ForwardDiff.gradient(vol, x0) ≈ Zygote.gradient(vol, x0)[1] rtol = 1e-12

    # frame + semi-rigid joints: Dual sections (promotion) + Dual EndConditions
    mat = Material(200.0, 77.0, 1.0, 0.3)
    fsec = Section(mat, 1e4, 8e7, 3e7, 5e6)
    fn1 = Node([0.0, 0.0, 0.0], :fixed)
    fn2 = Node([3000.0, 0.0, 0.0], :free)
    fn3 = Node([6000.0, 0.0, 0.0], :pinned)
    b1 = FrameElement(fn1, fn2, fsec,
        EndConditions(EndSprings(Inf, Inf, 1e8, 1e8), rigid_end()), :b1; rollangle=0.0)
    b2 = FrameElement(fn2, fn3, fsec, :b2; rollangle=0.0)
    fmodel = Asap.Model([fn1, fn2, fn3], AbstractElement{Float64}[b1, b2],
        AbstractLoad{Float64}[NodeForce(fn2, [0.0, -100.0, 0.0])])
    fp = OptParams(fmodel, AbstractVariable[
        JointVariable(b1, :start, 1e8, 1e5, 1e12),
        AreaVariable(b2, 1e4, 1e2, 1e5),
        SpatialVariable(fn2, 0.0, -500.0, 500.0, :Y)])
    fobj(x) = compliance(solve_structure(x, fp), fp)
    xf = copy(fp.values)
    @test ForwardDiff.gradient(fobj, xf) ≈ Zygote.gradient(fobj, xf)[1] rtol = 1e-10

    # network (FDM) path: Dual force densities through the multi-RHS solve
    ns = [Asap.FDMnode([Float64(i), Float64(j), 0.0], !(i in (0, 2) && j in (0, 2)))
          for j in 0:2 for i in 0:2]
    idx(i, j) = 3j + i + 1
    els = Asap.FDMelement[]
    for j in 0:2, i in 0:1
        push!(els, Asap.FDMelement(ns, idx(i, j), idx(i + 1, j), 1.0))
    end
    for j in 0:1, i in 0:2
        push!(els, Asap.FDMelement(ns, idx(i, j), idx(i, j + 1), 1.0))
    end
    loads = [Asap.FDMload(n, [0.0, 0.0, -1.0]) for n in ns if all(n.fixity)]
    network = Asap.Network(ns, els, loads)
    qv = QVariable(els[1], 1.5, 0.1, 10.0)
    np = NetworkOptParams(network, AbstractVariable[qv,
        CoupledVariable(els[2], qv), QVariable(els[7], 2.0, 0.1, 10.0)])
    nobj(x) = begin
        r = solve_network(x, np)
        sum(abs2, member_forces(r)) / 100
    end
    xn = copy(np.values)
    @test ForwardDiff.gradient(nobj, xn) ≈ Zygote.gradient(nobj, xn)[1] rtol = 1e-10
end

@testset "solver backends through OptParams" begin
    using LinearSolve

    model, vars = testbed()
    p0 = OptParams(model, vars)
    x = copy(p0.values) .+ [0.3, -0.1, 1e-3, 2e-3]
    u0 = solve_structure(x, p0).U
    obj(pp) = x -> compliance(solve_structure(x, pp), pp)
    g0 = Zygote.gradient(obj(p0), x)[1]

    for solver in (KLUFactorization(), KrylovJL_CG(), Asap.CachedSolver())
        p = OptParams(model, vars; solver = solver)
        # values, reverse gradients, forward gradients: all backend-invariant
        @test solve_structure(x, p).U ≈ u0 rtol = 1e-8
        @test Zygote.gradient(obj(p), x)[1] ≈ g0 rtol = 1e-6
        @test ForwardDiff.gradient(obj(p), x) ≈ g0 rtol = 1e-6
    end

    # CachedSolver shares ONE factorization across an entire ForwardDiff
    # Jacobian (all chunks see the same design values)
    cs = Asap.CachedSolver()
    pc = OptParams(model, vars; solver = cs)
    cstr(y) = axial_stress(solve_structure(y, pc), pc)
    # chunk size 2 over 4 vars → 2 chunks; the second must reuse the
    # factorization (same design values)
    cfg = ForwardDiff.JacobianConfig(cstr, x, ForwardDiff.Chunk{2}())
    J = ForwardDiff.jacobian(cstr, x, cfg)
    @test J ≈ Zygote.jacobian(y -> axial_stress(solve_structure(y, p0), p0), x)[1] rtol = 1e-8
    @test cs.hits > 0                       # chunks after the first reused it

    # network path accepts a solver too
    ns = [Asap.FDMnode([Float64(i), Float64(j), 0.0], !(i in (0, 2) && j in (0, 2)))
          for j in 0:2 for i in 0:2]
    idx(i, j) = 3j + i + 1
    els = Asap.FDMelement[]
    for j in 0:2, i in 0:1
        push!(els, Asap.FDMelement(ns, idx(i, j), idx(i + 1, j), 1.0))
    end
    for j in 0:1, i in 0:2
        push!(els, Asap.FDMelement(ns, idx(i, j), idx(i, j + 1), 1.0))
    end
    loads = [Asap.FDMload(n, [0.0, 0.0, -1.0]) for n in ns if all(n.fixity)]
    network = Asap.Network(ns, els, loads)
    qv = QVariable(els[1], 1.5, 0.1, 10.0)
    np0 = NetworkOptParams(network, AbstractVariable[qv, QVariable(els[7], 2.0, 0.1, 10.0)])
    npc = NetworkOptParams(network, AbstractVariable[qv, QVariable(els[7], 2.0, 0.1, 10.0)];
        solver = Asap.CachedSolver())
    xq = copy(np0.values)
    r0 = solve_network(xq, np0)
    rc = solve_network(xq, npc)
    @test [rc.X rc.Y rc.Z] ≈ [r0.X r0.Y r0.Z] rtol = 1e-10
    nobj(pp) = y -> sum(abs2, member_forces(solve_network(y, pp))) / 100
    @test Zygote.gradient(nobj(npc), xq)[1] ≈ Zygote.gradient(nobj(np0), xq)[1] rtol = 1e-8
end

@testset "kitchen sink: ALL variable types mixed in one problem" begin
    # frame + truss model carrying every current variable type at once:
    # spatial (X and Z), truss areas (grouped), frame area, joint stiffness
    # (with a mirrored tuple coupling) — gradients must agree across
    # Zygote, ForwardDiff, and finite differences
    mat = Material(200.0, 77.0, 1.0, 0.3)
    fsec = Section(mat, 1e4, 8e7, 3e7, 5e6)
    tsec = Section(mat, 1e3)

    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([3000.0, 0.0, 0.0], :free)
    n3 = Node([6000.0, 0.0, 0.0], :free)
    n4 = Node([9000.0, 0.0, 0.0], :pinned)
    n5 = Node([3000.0, 0.0, 3000.0], :pinned)
    n6 = Node([6000.0, 0.0, 3000.0], :pinned)

    b1 = FrameElement(n1, n2, fsec, EndConditions(EndSprings(Inf, Inf, 1e8, 1e8), rigid_end()), :b1; rollangle=0.0)
    b2 = FrameElement(n2, n3, fsec, :b2; rollangle=0.0)
    b3 = FrameElement(n3, n4, fsec, EndConditions(rigid_end(), EndSprings(Inf, Inf, 1e8, 1e8)), :b3; rollangle=0.0)
    t1 = TrussElement(n2, n5, tsec, :tie)
    t2 = TrussElement(n3, n6, tsec, :tie)

    model = Asap.Model([n1, n2, n3, n4, n5, n6],
        AbstractElement{Float64}[b1, b2, b3, t1, t2],
        AbstractLoad{Float64}[NodeForce(n2, [0.0, -100.0, 0.0]),
            NodeForce(n3, [0.0, -80.0, 20.0])])

    sv = SpatialVariable(n5, 0.0, -500.0, 500.0, :Z)
    av_t = AreaVariable(t1, 1e3, 1e2, 1e4)
    jv = JointVariable(b1, :start, 1e8, 1e5, 1e12)
    vars = AbstractVariable[
        sv,
        CoupledVariable(n6, sv, -1.0),            # mirrored spatial partner
        SpatialVariable(n2, 0.0, -400.0, 400.0, :X),
        av_t,
        CoupledVariable(t2, av_t),                # grouped truss areas
        AreaVariable(b2, 1e4, 1e3, 1e5),          # FRAME area (Dual section promotion)
        jv,
        CoupledVariable((b3, :end), jv, 1.0),     # mirrored joint pair
    ]
    p = OptParams(model, vars)
    @test length(p.values) == 5                   # independents only

    x = copy(p.values) .+ [50.0, -30.0, 2e2, 3e3, 1e8]
    obj(y) = begin
        res = solve_structure(y, p)
        compliance(res, p) + 1e-8 * sum(abs2, axial_stress(res, p))
    end
    g_zy = Zygote.gradient(obj, x)[1]
    g_fw = ForwardDiff.gradient(obj, x)
    g_fd = FiniteDifferences.grad(FDM, obj, x)[1]
    @test g_zy ≈ g_fw rtol = 1e-9
    # FD truncation at the 1e8 joint-stiffness scale dominates the tolerance
    @test g_zy ≈ g_fd rtol = 1e-3
    @test all(abs.(g_zy) .> 0)                    # every variable type is live

    # materialized model agrees with the pure path at this mixed design
    m2 = updatemodel(p, x)
    @test solve_structure(x, p).U ≈ m2.results.u rtol = 1e-10
end

@testset "implicit differentiation (solution_tangents)" begin
    # mixed spatial + grouped areas (testbed truss): every Jacobian must
    # match full-pipeline ForwardDiff AND Zygote to machine precision
    model, vars = testbed()
    p = OptParams(model, vars)
    x = copy(p.values) .+ [0.3, -0.15, 2e-3, 4e-3]

    t = solution_tangents(x, p)
    @test t.res.U ≈ solve_structure(x, p).U rtol = 1e-12

    dU_fw = ForwardDiff.jacobian(y -> solve_structure(y, p).U, x)
    @test t.dU ≈ dU_fw rtol = 1e-10

    Jσ = axial_stress_jacobian(t, p)
    Jσ_fw = ForwardDiff.jacobian(y -> axial_stress(solve_structure(y, p), p), x)
    Jσ_zy = Zygote.jacobian(y -> axial_stress(solve_structure(y, p), p), x)[1]
    @test Jσ ≈ Jσ_fw rtol = 1e-10
    @test Jσ ≈ Jσ_zy rtol = 1e-10

    JN = axial_force_jacobian(t, p)
    JN_fw = ForwardDiff.jacobian(y -> axial_force(solve_structure(y, p), p), x)
    @test JN ≈ JN_fw rtol = 1e-10

    # kitchen-sink frame: spatial + truss/frame areas + joints + couplings
    mat = Material(200.0, 77.0, 1.0, 0.3)
    fsec = Section(mat, 1e4, 8e7, 3e7, 5e6)
    tsec = Section(mat, 1e3)
    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([3000.0, 0.0, 0.0], :free)
    n3 = Node([6000.0, 0.0, 0.0], :free)
    n4 = Node([9000.0, 0.0, 0.0], :pinned)
    n5 = Node([3000.0, 0.0, 3000.0], :pinned)
    n6 = Node([6000.0, 0.0, 3000.0], :pinned)
    b1 = FrameElement(n1, n2, fsec, EndConditions(EndSprings(Inf, Inf, 1e8, 1e8), rigid_end()), :b1; rollangle=0.0)
    b2 = FrameElement(n2, n3, fsec, :b2; rollangle=0.0)
    b3 = FrameElement(n3, n4, fsec, EndConditions(rigid_end(), EndSprings(Inf, Inf, 1e8, 1e8)), :b3; rollangle=0.0)
    t1 = TrussElement(n2, n5, tsec, :tie)
    t2 = TrussElement(n3, n6, tsec, :tie)
    fmodel = Asap.Model([n1, n2, n3, n4, n5, n6],
        AbstractElement{Float64}[b1, b2, b3, t1, t2],
        AbstractLoad{Float64}[NodeForce(n2, [0.0, -100.0, 0.0]),
            NodeForce(n3, [0.0, -80.0, 20.0])])
    sv = SpatialVariable(n5, 0.0, -500.0, 500.0, :Z)
    av_t = AreaVariable(t1, 1e3, 1e2, 1e4)
    jv = JointVariable(b1, :start, 1e8, 1e5, 1e12)
    fp = OptParams(fmodel, AbstractVariable[
        sv, CoupledVariable(n6, sv, -1.0),
        SpatialVariable(n2, 0.0, -400.0, 400.0, :X),
        av_t, CoupledVariable(t2, av_t),
        AreaVariable(b2, 1e4, 1e3, 1e5),
        jv, CoupledVariable((b3, :end), jv, 1.0)])
    xf = copy(fp.values) .+ [50.0, -30.0, 2e2, 3e3, 1e8]

    tf = solution_tangents(xf, fp)
    dU_fw = ForwardDiff.jacobian(y -> solve_structure(y, fp).U, xf)
    @test tf.dU ≈ dU_fw rtol = 1e-9

    Jσf = axial_stress_jacobian(tf, fp)
    Jσf_fw = ForwardDiff.jacobian(y -> axial_stress(solve_structure(y, fp), fp), xf)
    @test Jσf ≈ Jσf_fw rtol = 1e-9

    # CachedSolver path shares the factorization
    cs = Asap.CachedSolver()
    pc = OptParams(model, vars; solver = cs)
    tc = solution_tangents(x, pc)
    @test tc.dU ≈ t.dU rtol = 1e-10
end

@testset "directional rail SpatialVariables" begin
    mat = Material(200e6, 1.0, 80.0, 0.3)
    sec = Section(mat, 1e-2)
    rot = [true, true, true]
    n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = Node([4.0, 0.0, 0.0], vcat([true, true, false], rot))
    n3 = Node([8.0, 0.0, 0.0], vcat([false, false, false], rot))
    n4 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot), :top)
    n5 = Node([6.0, 3.0, 0.0], vcat([true, true, false], rot), :top)
    els = AbstractElement{Float64}[
        TrussElement(n1, n2, sec), TrussElement(n2, n3, sec), TrussElement(n4, n5, sec),
        TrussElement(n1, n4, sec), TrussElement(n4, n2, sec),
        TrussElement(n2, n5, sec), TrussElement(n5, n3, sec)]
    model = Model([n1, n2, n3, n4, n5], els,
        AbstractLoad{Float64}[NodeForce(n2, [0.0, -50.0, 0.0])])

    # rail along an UNNORMALIZED diagonal (constructor must normalize),
    # with a mirrored coupled partner and a mixed-in area variable
    rail = SpatialVariable(n4, [2.0, 2.0, 0.0], 0.0, -1.0, 1.5)
    @test rail.direction ≈ [1 / √2, 1 / √2, 0.0]
    @test SpatialVariable(n4, [0.0, 3.0, 0.0], -1.0, 1.0).value == 0.0  # legacy 4-arg form

    av = AreaVariable(els[1], 1e-2, 1e-4, 5e-2)
    p = OptParams(model, AbstractVariable[
        rail, CoupledVariable(n5, rail, -1.0), av])
    @test length(p.values) == 2                    # rail + area (coupling shares)
    x = copy(p.values) .+ [0.4, 3e-3]

    # materialization: node moved by value·unit(direction), mirror negated
    m2 = updatemodel(p, x)
    @test m2.nodes[4].position ≈ [2.0, 3.0, 0.0] .+ 0.4 .* [1 / √2, 1 / √2, 0.0]
    @test m2.nodes[5].position ≈ [6.0, 3.0, 0.0] .- 0.4 .* [1 / √2, 1 / √2, 0.0]
    @test solve_structure(x, p).U ≈ m2.results.u rtol = 1e-10

    # gradients through every path
    obj(y) = compliance(solve_structure(y, p), p)
    g_zy = Zygote.gradient(obj, x)[1]
    @test g_zy ≈ ForwardDiff.gradient(obj, x) rtol = 1e-10
    @test g_zy ≈ FiniteDifferences.grad(FDM, obj, x)[1] rtol = 1e-6

    # implicit differentiation inherits rails through the scatter maps
    t = solution_tangents(x, p)
    @test t.dU ≈ ForwardDiff.jacobian(y -> solve_structure(y, p).U, x) rtol = 1e-10
    @test axial_stress_jacobian(t, p) ≈
          Zygote.jacobian(y -> axial_stress(solve_structure(y, p), p), x)[1] rtol = 1e-10

    # axis constructors are the unit-vector special case (incl. lowercase)
    @test SpatialVariable(n4, 0.0, -1.0, 1.0, :Z).direction == [0.0, 0.0, 1.0]
    @test SpatialVariable(n4, 0.0, -1.0, 1.0, :y).direction == [0.0, 1.0, 0.0]
    @test_throws ArgumentError SpatialVariable(n4, 0.0, -1.0, 1.0, :W)
end
