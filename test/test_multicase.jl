# Multi-case solves (LoadCases / solve_cases) and mirrored / directional couplings.
# Motivating use: code-checked truss optimization over many load cases and
# combinations (CimbraWoodCore's TrussOptProblem built one OptParams per case).

# the design-dependent-loads portal, loads tagged by case; `keep` filters them
function mc_portal(keep=(:dead, :live, :wind))
    mat = Material(200e6, 80e6, 8.0, 0.3)
    col = Section(mat, 1.2e-2, 3e-4, 1e-4, 2e-6)
    bm = Section(mat, 8e-3, 2.5e-4, 5e-5, 1e-6)
    planar = [true, true, false, false, false, true]
    n1 = Node([0.0, 0.0, 0.0], :fixed); n2 = Node([0.0, 3.0, 0.0], planar)
    n3 = Node([6.0, 4.0, 0.0], planar); n4 = Node([6.0, 0.0, 0.0], :fixed)
    sp = EndSprings(Inf, Inf, 5e4, 5e4)
    c1 = FrameElement(n1, n2, col, :col; rollangle=0.0)
    b = FrameElement(n2, n3, bm, EndConditions(sp, sp), :beam; rollangle=0.0)
    c2 = FrameElement(n4, n3, col, :col; rollangle=0.0)
    all_loads = AbstractLoad{Float64}[
        LineLoad(b, [0.0, -8.0, 0.0]; case=:dead),
        SelfWeight(c1; g=[0.0, -9.81, 0.0], case=:dead),
        SelfWeight(c2; g=[0.0, -9.81, 0.0], case=:dead),
        PointLoad(b, 0.3, [0.0, -5.0, 0.0]; case=:live),
        LineLoad(b, [0.0, -3.0, 0.0]; case=:live),
        NodeForce(n2, [4.0, 0.0, 0.0]; case=:wind),
        NodeForce(n3, [2.5, 0.0, 0.0]; case=:wind)]
    loads = AbstractLoad{Float64}[l for l in all_loads if l.case in keep]
    model = Asap.Model([n1, n2, n3, n4], AbstractElement{Float64}[c1, b, c2], loads)
    y2 = SpatialVariable(n2, 0.0, -1.0, 1.0, :Y)
    x3 = SpatialVariable(n3, 0.0, -1.0, 1.0, :X)
    Ac = AreaVariable(c1, 1.2e-2, 1e-3, 5e-2)
    Ib = SectionVariable(b, 2.5e-4, 1e-5, 1e-3, :Ix)
    jv = JointVariable(b, :start, 5e4, 1e2, 1e8)
    vars = AbstractVariable[y2, x3, Ac, CoupledVariable(c2, Ac), Ib, jv, CoupledVariable((b, :end), jv)]
    return model, vars
end

@testset "multi-case solves (LoadCases, solve_cases)" begin
    for loads in (:design, :frozen)
        p = OptParams(mc_portal()...; loads)
        x0 = copy(p.values)
        x1 = x0 .* [1, 1, 1.3, 0.8, 1.5] .+ [0.4, -0.3, 0, 0, 0]
        r_all = solve_structure(x1, p)
        lc = LoadCases(p)
        @test lc.names == [:dead, :live, :wind]
        @test solve_structure(x1, p).U == r_all.U          # constructing LoadCases restored the cache
        rc = solve_cases(x1, p, lc)
        @test size(rc.U) == (length(r_all.U), 3)
        # superposition: the cases add up to the all-loads solve
        @test vec(sum(rc.U; dims=2)) ≈ r_all.U rtol = 1e-10
        @test vec(sum(rc.F; dims=2)) ≈ r_all.F rtol = 1e-10
        # parity with one OptParams per case (the old workaround)
        Nmc = axial_force(rc, p, lc)
        cmc = compliance(rc, p, lc)
        for (k, c) in enumerate(lc.names)
            pk = OptParams(mc_portal((c,))...; loads)
            rk = solve_structure(x1, pk)
            @test rc.U[:, k] ≈ rk.U rtol = 1e-10
            @test rc.F[:, k] ≈ rk.F rtol = 1e-10
            @test cmc[k] ≈ compliance(rk, pk) rtol = 1e-10
            @test Nmc[:, k] ≈ axial_force(rk, pk) rtol = 1e-9 atol = 1e-12
            @test axial_stress(rc, p, lc)[:, k] ≈ axial_stress(rk, pk) rtol = 1e-9 atol = 1e-9
        end
        # gradients: Zygote = ForwardDiff = Σ of the verified single-case gradients
        w = [1.0, 0.7, 0.3]
        obj(x) = sum(w .* compliance(solve_cases(x, p, lc), p, lc))
        g_zy = Zygote.gradient(obj, x0)[1]
        g_fd = ForwardDiff.gradient(obj, x0)
        g_ref = sum(w[k] * ForwardDiff.gradient(x -> (pk = OptParams(mc_portal((c,))...; loads);
                                                       compliance(solve_structure(x, pk), pk)), x0)
                    for (k, c) in enumerate(lc.names))
        @test g_zy ≈ g_fd rtol = 1e-9
        @test g_fd ≈ g_ref rtol = 1e-9
        # constraint Jacobian of the per-case axial forces (forward mode)
        JN = ForwardDiff.jacobian(x -> vec(axial_force(solve_cases(x, p, lc), p, lc)), x1)
        pw = OptParams(mc_portal((:live,))...; loads)
        @test JN[4:6, :] ≈ ForwardDiff.jacobian(x -> axial_force(solve_structure(x, pw), pw), x1) rtol = 1e-8 atol = 1e-12
        # a combination is a linear map of the case columns
        @test rc.U * [1.2, 1.6, 0.0] ≈ 1.2 .* rc.U[:, 1] .+ 1.6 .* rc.U[:, 2]
    end
    # selection, errors, show
    p = OptParams(mc_portal()...)
    lc2 = LoadCases(p; cases=[:wind, :dead])
    @test lc2.names == [:wind, :dead] && size(solve_cases(p.values, p, lc2).U, 2) == 2
    @test_throws ArgumentError LoadCases(p; cases=[:snow])
    @test_throws ArgumentError LoadCases(p; cases=[:dead, :dead])
    @test occursin(":wind, :dead", sprint(show, MIME("text/plain"), lc2))
end

@testset "mirrored and directional couplings" begin
    mat = Material(200e6, 1.0, 80.0, 0.3)
    sec = Section(mat, 1e-2)
    rot = [true, true, true]
    n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = Node([4.0, 0.0, 0.0], vcat([true, true, false], rot))
    n3 = Node([8.0, 0.0, 0.0], vcat([true, false, false], rot))
    n4 = Node([2.0, 1.5, 0.0], vcat([true, true, false], rot))
    n5 = Node([6.0, 1.5, 0.0], vcat([true, true, false], rot))
    n6 = Node([4.0, 3.0, 0.0], vcat([true, true, false], rot))
    els = AbstractElement{Float64}[
        TrussElement(n1, n2, sec), TrussElement(n2, n3, sec),
        TrussElement(n1, n4, sec), TrussElement(n4, n6, sec), TrussElement(n6, n5, sec), TrussElement(n5, n3, sec),
        TrussElement(n4, n2, sec), TrussElement(n2, n6, sec), TrussElement(n2, n5, sec)]
    loads = AbstractLoad{Float64}[NodeForce(n4, [0.0, -20.0, 0.0]), NodeForce(n5, [0.0, -20.0, 0.0]),
        NodeForce(n6, [0.0, -30.0, 0.0])]
    model = Model([n1, n2, n3, n4, n5, n6], els, loads)

    # n4 slides along the pitched top chord (rail toward the ridge); n5 mirrors it across x = 4
    rail = SpatialVariable(n4, [2.0, 1.5, 0.0], 0.0, -1.0, 1.0)
    mir = CoupledVariable(n5, rail; mirror=:X)
    @test mir.direction ≈ [-0.8, 0.6, 0.0] && mir.factor == 1.0
    @test CoupledVariable(n5, rail; mirror=[2.0, 0.0, 0.0]).direction ≈ mir.direction
    p = OptParams(model, AbstractVariable[rail, mir])
    @test length(p.values) == 1
    x = [0.37]
    m2 = updatemodel(p, x)
    a, b = m2.nodes[4].position, m2.nodes[5].position
    @test b[1] ≈ 8.0 - a[1] && b[2] ≈ a[2] && a ≈ [2.0, 1.5, 0.0] .+ 0.37 .* [0.8, 0.6, 0.0]
    # symmetric structure, symmetric load → symmetric response (n1 pin, n3 roller)
    r = solve_structure(x, p)
    u(i) = r.U[6(i-1)+1:6(i-1)+2]
    @test u(5)[2] ≈ u(4)[2] rtol = 1e-10
    @test r.U ≈ m2.results.u rtol = 1e-10
    obj(y) = compliance(solve_structure(y, p), p)
    g = Zygote.gradient(obj, x)[1]
    @test g ≈ ForwardDiff.gradient(obj, x) rtol = 1e-10
    @test g ≈ FiniteDifferences.grad(FDM, obj, x)[1] rtol = 1e-6
    t = solution_tangents(x, p)
    @test t.dU ≈ ForwardDiff.jacobian(y -> solve_structure(y, p).U, x) rtol = 1e-10

    # explicit target direction (normalized), factor applied along it
    dv = CoupledVariable(n6, rail, 0.5; direction=[0.0, 2.0, 0.0])
    @test dv.direction == [0.0, 1.0, 0.0]
    p2 = OptParams(model, AbstractVariable[rail, dv])
    @test updatemodel(p2, x).nodes[6].position ≈ [4.0, 3.0 + 0.5 * 0.37, 0.0]
    # the default coupling is unchanged (parent's direction)
    @test CoupledVariable(n5, rail).direction === nothing
    @test occursin("own direction", sprint(show, MIME("text/plain"), mir))
    # guards
    av = AreaVariable(els[1], 1e-2, 1e-4, 5e-2)
    @test_throws ArgumentError CoupledVariable(els[2], av; mirror=:X)
    @test_throws ArgumentError CoupledVariable(n5, rail; mirror=:X, direction=[1.0, 0.0, 0.0])
    @test_throws ArgumentError CoupledVariable(n5, rail; direction=[0.0, 0.0, 0.0])
    @test_throws ArgumentError CoupledVariable(n5, rail; mirror=:W)
end
