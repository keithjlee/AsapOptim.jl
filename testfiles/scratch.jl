"""
    displacement(values::Vector{Float64}, p::TrussOptParams)
    
Get the vector of DOF displacements given a set of design variables and problem parameters. This function is the basis of ALL subsequent structural analysis
"""
function f1(values::Vector{Float64}, p::TrussOptParams)
    
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
    rotmats = AsapOptim.getRmatrices(elementvecsnormalized, p)

    # kₑ
    klocs = AsapOptim.klocal.(p.E, Anew, elementlengths)

    # Kₑ
    kglobs = AsapOptim.getglobalks(rotmats, klocs, p)

    # K
    K = AsapOptim.assembleglobalK(kglobs, p)

    # K⁻¹P
    disp = AsapOptim.solveU(K, p)

    # U
    U = AsapOptim.replacevalues(zeros(p.n), p.freeids, disp)

    F = getindex.(rotmats .* kglobs .* [U[id] for id in p.dofids], 2)

    abs.(F)' * elementlengths
end

function f2(values::Vector{Float64}, p::TrussOptParams)
    
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
    rotmats = AsapOptim.getRmatrices(elementvecsnormalized, p)

    # kₑ
    klocs = AsapOptim.klocal.(p.E, Anew, elementlengths)

    # Kₑ
    kglobs = AsapOptim.getglobalks(rotmats, klocs, p)

    # K
    K = AsapOptim.assembleglobalK(kglobs, p)

    # K⁻¹P
    disp = AsapOptim.solveU(K, p)

    # U
    U = AsapOptim.replacevalues(zeros(p.n), p.freeids, disp)

    F = AsapOptim.axialforces(U, kglobs, rotmats, p)

    abs.(F)' * elementlengths
end

@time f1(vals, problem)
@time f2(vals, problem)

@time g1 = Zygote.gradient(var -> f1(var, problem), vals)[1];
@time g2 = Zygote.gradient(var -> f2(var, problem), vals)[1];

norm(g1 .- g2)

gdiff = g1 .- g2
inz = findall(abs.(gdiff) .>= 1e-5)

varnz = vars[inz]