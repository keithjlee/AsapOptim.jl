#baseline
function obj2(values::Vector{Float64}, p::TrussOptParams)
    indexer = p.indexer

    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])
    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])

    ks = [kglobal(Xnew, Ynew, Znew, e, a, id) for (e, a, id) in zip(p.E, Anew, p.nodeids)]

    K = assembleglobalK(ks, p)

    U = solveU(K, p)

    return U' * p.P[p.freeids]
end

@time o2 = obj2(vals, problem);
@time g2 = Zygote.gradient(var -> obj2(var, problem), vals)[1];

#optim
function obj_test(values::Vector{Float64}, p::TrussOptParams)
    #populate values
    indexer = p.indexer
    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])
    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])

    # vₑ
    # elementvecs = AsapOptim.getevecs(Xnew, Ynew, Znew, p)
    elementvecs = p.C * [Xnew Ynew Znew]

    # Lₑ
    # elementlengths = AsapOptim.getlengths(elementvecs)
    elementlengths = norm.(eachrow(elementvecs))

    # vnₑ
    # elementvecsnormalized = AsapOptim.getnormalizedevecs(elementvecs, elementlengths)
    elementvecsnormalized = elementvecs ./ elementlengths

    # Γ
    rotmats = AsapOptim.Rtruss.(eachrow(elementvecsnormalized))

    # kₑ
    klocs = AsapOptim.klocal.(p.E, Anew, elementlengths)

    # Kₑ
    # kglobs = AsapOptim.getglobalks(rotmats, klocs)
    kglobs = transpose.(rotmats) .* klocs .* rotmats

    # K
    K = AsapOptim.assembleglobalK(kglobs, p)

    # K⁻¹P
    disp = AsapOptim.solveU(K, p)

    disp' * p.P[p.freeids]
end

@time otest = obj_test(vals, problem);
@time gtest = Zygote.gradient(var -> obj_test(var, problem), vals)[1];

@show  o2 - otest
@show norm(g2 .- gtest)