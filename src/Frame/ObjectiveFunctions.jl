"""
    displacement(values::Vector{Float64}, p::FrameOptParams)
    
Get the vector of DOF displacements given a set of design variables and problem parameters. This function is the basis of ALL subsequent structural analysis
"""
function displacement(values::Vector{Float64}, p::FrameOptParams)
    
    #populate values
    Xnew = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Ynew = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Znew = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    Anew = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])
    Iznew = replacevalues(p.Iz, p.indexer.iIz, values[p.indexer.iIzg])
    Iynew = replacevalues(p.Iy, p.indexer.iIy, values[p.indexer.iIyg])

    # vₑ
    elementvecs = getevecs(Xnew, Ynew, Znew, p)

    # Lₑ
    elementlengths = getlengths(elementvecs, p)

    # vnₑ
    elementvecsnormalized = getnormalizedevecs(elementvecs, elementlengths)

    # Γ
    rotmats = getRmatrices(elementvecsnormalized, p)

    # kₑ
    klocs = ktruss.(p.E, Anew, elementlengths)

    # Kₑ
    kglobs = getglobalks(rotmats, klocs, p)

    # K
    K = assembleglobalK(kglobs, p)

    # K⁻¹P
    disp = solveU(K, p)

    # U
    replacevalues(zeros(p.n), p.freeids, disp)
end