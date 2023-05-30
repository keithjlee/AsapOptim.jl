function obj_test(values::Vector{Float64}, p::TrussOptParams)
    #populate values
    Xnew = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Ynew = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Znew = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    Anew = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    elementvecs = AsapOptim.getevecs(Xnew, Ynew, Znew, p)

    # Lₑ
    elementlengths = AsapOptim.getlengths(elementvecs, p)

    # vnₑ
    elementvecsnormalized = AsapOptim.getnormalizedevecs(elementvecs, elementlengths)

    # Γ
    # rotmats = AsapOptim.Rtruss.(eachrow(elementvecsnormalized))
    rotmats = AsapOptim.getRmatrices(elementvecsnormalized, p)

    # kₑ
    # klocs = AsapOptim.getlocalks(p.E, Anew, elementlengths)
    klocs = AsapOptim.klocal.(p.E, Anew, elementlengths)

    # Kₑ
    kglobs = AsapOptim.getglobalks(rotmats, klocs, p)

    # K
    K = AsapOptim.assembleglobalK(kglobs, p)

    # K⁻¹P
    disp = solveU(K, p)

    # U
    U = replacevalues(zeros(p.n), p.freeids, disp)

    U' * p.P
end;

@time otest = obj_test(vals, problem);
@time gtest = Zygote.gradient(var -> obj_test(var, problem), vals)[1];

@show norm(gtest .- g1)


m1 = rand(2,6); m2 = rand(2,6); m3 = rand(2,6)
ts = [rand(2,6) for _ = 1:3000]

#
@time r1 = [[dot(t, m1) for t in ts] [dot(t, m2) for t in ts] [dot(t, m3) for t in ts]];

@time begin
    r2 = zeros(length(ts), 3)

    @inbounds for i in axes(ts, 1)
        r2[i, 1] = dot(ts[i], m1)
        r2[i, 2] = dot(ts[i], m2)
        r2[i, 3] = dot(ts[i], m3)
    end
end;

tt = rand(1000, 3);

@time tt1 = zero(tt);
@time tt2 = zeros(size(tt, 1), 3);