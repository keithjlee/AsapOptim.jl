function F1(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = AsapOptim.addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Y = AsapOptim.addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Z = AsapOptim.addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    A = AsapOptim.replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    v = AsapOptim.getevecs(X, Y, Z, p)

    # Lₑ
    l = AsapOptim.getlengths(v)

    # vnₑ
    n = AsapOptim.getnormalizedevecs(v, l)

    # Γ
    Γ = AsapOptim.getRmatrices(n)

    # kₑ
    kₑ = AsapOptim.ktruss.(p.E, A, l)

    # Kₑ
    Kₑ = AsapOptim.getglobalks(Γ, kₑ)

    # K
    K = AsapOptim.assembleglobalK(Kₑ, p)

    # K⁻¹P
    u = AsapOptim.solveU(K, p)

    # U
    U = AsapOptim.replacevalues(zeros(p.n), p.freeids, u)

    F = AsapOptim.Faxial(U, Kₑ, Γ, p)

    maximum(F) - minimum(F)
end

function F2(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = AsapOptim.addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Y = AsapOptim.addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Z = AsapOptim.addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    A = AsapOptim.replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    v = AsapOptim.getevecs(X, Y, Z, p)

    # Lₑ
    l = AsapOptim.getlengths(v)

    # vnₑ
    n = AsapOptim.getnormalizedevecs(v, l)

    # Γ
    Γ = AsapOptim.getRmatrices(n)

    # kₑ
    kₑ = AsapOptim.ktruss.(p.E, A, l)

    # Kₑ
    Kₑ = AsapOptim.getglobalks(Γ, kₑ)

    # K
    K = AsapOptim.assembleglobalK(Kₑ, p)

    # K⁻¹P
    u = AsapOptim.solveU(K, p)

    # U
    U = AsapOptim.replacevalues(zeros(p.n), p.freeids, u)

    # Store values for continuity in gradients
    tr = AsapOptim.TrussResults(X,
        Y,
        Z,
        A,
        l,
        Kₑ,
        Γ,
        U)

    F = AsapOptim.axialforce(tr, p)

    maximum(F) - minimum(F)
end

function F3(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = AsapOptim.addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Y = AsapOptim.addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Z = AsapOptim.addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    A = AsapOptim.replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    v = AsapOptim.getevecs(X, Y, Z, p)

    # Lₑ
    l = AsapOptim.getlengths(v)

    # vnₑ
    n = AsapOptim.getnormalizedevecs(v, l)

    # Γ
    Γ = AsapOptim.getRmatrices(n)

    # kₑ
    kₑ = AsapOptim.ktruss.(p.E, A, l)

    # Kₑ
    Kₑ = AsapOptim.getglobalks(Γ, kₑ)

    # K
    K = AsapOptim.assembleglobalK(Kₑ, p)

    # K⁻¹P
    u = AsapOptim.solveU(K, p)

    # U
    U = AsapOptim.replacevalues(zeros(p.n), p.freeids, u)

    F = AsapOptim.Faxial2(U, Kₑ, Γ, p)

    maximum(F) - minimum(F)
end


function Fbaseline(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = AsapOptim.addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Y = AsapOptim.addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Z = AsapOptim.addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    A = AsapOptim.replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    v = AsapOptim.getevecs(X, Y, Z, p)

    # Lₑ
    l = AsapOptim.getlengths(v)

    # vnₑ
    n = AsapOptim.getnormalizedevecs(v, l)

    # Γ
    Γ = AsapOptim.getRmatrices(n)

    # kₑ
    kₑ = AsapOptim.ktruss.(p.E, A, l)

    # Kₑ
    Kₑ = AsapOptim.getglobalks(Γ, kₑ)

    # K
    K = AsapOptim.assembleglobalK(Kₑ, p)

    # K⁻¹P
    u = AsapOptim.solveU(K, p)

    # U
    U = AsapOptim.replacevalues(zeros(p.n), p.freeids, u)

    Ue = [U[id] for id in p.dofids]

    Fvecs = Γ .* Kₑ .* Ue

    F = getindex.(Fvecs, 2)

    maximum(F) - minimum(F)
end

@time o1 = F1(vals, problem)
@time o2 = F2(vals, problem)
@time o3 = F3(vals, problem)
@time o = Fbaseline(vals, problem)

@time g1 = Zygote.gradient(x -> F1(x, problem), vals)[1];
@time g2 = Zygote.gradient(x -> F2(x, problem), vals)[1];
@time g3 = Zygote.gradient(x -> F3(x, problem), vals)[1];
@time g = Zygote.gradient(x -> Fbaseline(x, problem), vals)[1];

@show norm(g1 .- g)
@show norm(g2 .- g)
@show norm(g3 .- g)

gdiff = g2 .- g3

vnz = vars[findall(gdiff .!= 0)]
getproperty.(vnz, :i)

u = rand(6)
k = rand(6,6)
r = rand(2,6)
f = rand(2)

kron(u', r)' * f