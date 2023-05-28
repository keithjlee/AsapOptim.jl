
#enzyme test

using Enzyme

function objEnzyme(u::Vector{Float64}, values::Vector{Float64}, p::TrussOptParams)

    indexer = p.indexer
    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])
    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])

    #element vectors
    vₑ = getevecs(Xnew, Ynew, Znew, p)

    #element lengths
    Lₑ = getlengths(vₑ)

    #normalized vecs
    nₑ = getnormalizedevecs(vₑ, Lₑ)

    #transformation matrices
    Γₑ = AsapOptim.Rtruss.(eachrow(nₑ))

    #local stiffness matrices
    kₑ = AsapOptim.klocal.(p.E, Anew, Lₑ)

    #global stiffness matrices
    Kₑ = getglobalks(Γₑ, kₑ)

    #global stiffness matrix
    K = AsapOptim.assembleglobalK(Kₑ, p)

    #solve displacement
    U = solveU(K, p)

    u[1] = U' * p.P[p.freeids]
end


x = vals
dx = zero(x)

y = [0.]
dy = [1.]

autodiff(Enzyme.Reverse, objEnzyme, Duplicated(y, dy), Duplicated(x, dx), Const(p))