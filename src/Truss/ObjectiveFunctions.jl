"""
    displacement(values::Vector{Float64}, p::TrussOptParams)
    
Get the vector of DOF displacements given a set of design variables and problem parameters. This function is the basis of ALL subsequent structural analysis
"""
function displacement(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    Xnew = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Ynew = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Znew = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    Anew = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    elementvecs = getevecs(Xnew, Ynew, Znew, p)

    # Lₑ
    elementlengths = getlengths(elementvecs)

    # vnₑ
    elementvecsnormalized = getnormalizedevecs(elementvecs, elementlengths)

    # Γ
    rotmats = getRmatrices(elementvecsnormalized)

    # kₑ
    klocs = ktruss.(p.E, Anew, elementlengths)

    # Kₑ
    kglobs = getglobalks(rotmats, klocs)

    # K
    K = assembleglobalK(kglobs, p)

    # K⁻¹P
    disp = solveU(K, p)

    # U
    replacevalues(zeros(p.n), p.freeids, disp)
end

"""
    compliance(u::Vector{Float64}, p::TrussOptParams)

Measure of strain energy for truss structures.
"""
function compliance(u::Vector{Float64}, p::TrussOptParams)
    u' * p.P
end