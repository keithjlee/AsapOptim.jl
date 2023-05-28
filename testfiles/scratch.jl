indexer = p.indexer

Xnew = addvalues(p.X, indexer.iX, vals[indexer.iXg])
Ynew = addvalues(p.Y, indexer.iY, vals[indexer.iYg])
Znew = addvalues(p.Z, indexer.iZ, vals[indexer.iZg])
Anew = replacevalues(p.A, indexer.iA, vals[indexer.iAg])

#element vectors
elementvecs = p.C * [Xnew Ynew Znew]

#element lengths
elementlengths = norm.(eachrow(elementvecs))

#normalized vecs
elementvecsnormalized = elementvecs ./ elementlengths

#Γ
@time rotmats = AsapOptim.Rtruss.(eachrow(elementvecsnormalized));
@time rotmats2 = getrs(elementvecsnormalized);

#kloc
klocs = AsapOptim.klocal.(p.E, Anew, elementlengths)

#kglob
kglobs = transpose.(rotmats) .* klocs .* rotmats

#global stiffness matrix
K = AsapOptim.assembleglobalK(kglobs, p)

#solve displacement
U = solveU(K, p)