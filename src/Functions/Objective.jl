"""
    solve(values::Vector{Float64}, p::TrussOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solvetruss(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX)
    Y = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY)
    Z = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ)
    A = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA)

    # vₑ
    v = getevecs(X, Y, Z, p)

    # Lₑ
    l = getlengths(v)

    # vnₑ
    n = getnormalizedevecs(v, l)

    # Γ
    Γ = Rtruss(n)

    # kₑ
    kₑ = ktruss.(p.E, A, l)

    # Kₑ
    Kₑ = getglobalks(Γ, kₑ)

    # K
    K = assembleglobalK(Kₑ, p)

    # K⁻¹P
    u = solveU(K, p)

    # U
    U = replacevalues(zeros(p.n), p.freeids, u)

    # Store values for continuity in gradients
    return TrussResults(X,
        Y,
        Z,
        A,
        l,
        Kₑ,
        Γ,
        U)
end

"""
    compliance(t::TrussResults, p::TrussOptParams)

Measure of strain energy for truss structures.
"""
function compliance(t::TrussResults, p::TrussOptParams)
    t.U' * p.P
end

"""
    variation(vals::Vector{Float64}; factor = 1.)
Penalize the distance between extrema of a set
"""
variation(vals::Vector{Float64}; factor = 1.) = -factor * reduce(-, extrema(vals))

"""
    maxpenalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
Penalize values above a threshold
"""
function maxpenalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
    Δ = vals .- threshold
    factor * sum(Δ .+ abs.(Δ))
end

"""
    minpenalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
Penalize values below a threshold
"""
function minpenalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
    Δ = threshold .- vals
    factor * sum(Δ .+ abs.(Δ))
end

"""
    volume(values::Vector{Float64}, p::TrussOptParams)
Extract the volume of a truss
"""
function volume(values::Vector{Float64}, p::TrussOptParams)

    #populate values
    X = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX)
    Y = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY)
    Z = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ)
    A = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA)

    # vₑ
    v = getevecs(X, Y, Z, p)

    # Lₑ
    l = getlengths(v)

    dot(A, l)
end

"""
    solvenetwork(values::Vector{Float64}, p::NetworkOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solvenetwork(values::Vector{Float64}, p::NetworkOptParams)
    
    #populate values
    X = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX)
    Y = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY)
    Z = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ)
    q = replacevalues(p.q, p.indexer.iQ, values[p.indexer.iQg] .* p.indexer.fQ)

    # fixed nodal positions
    xyz_f = [X[p.F] Y[p.F] Z[p.F]]

    # diagonal q matrix
    Q = diagm(q)

    #solve for free positions
    xyz_n = (p.Cn' * Q * p.Cn) \ (p.Pn - p.Cn' * Q * p.Cf * xyz_f)

    X2 = replacevalues(X, p.N, xyz_n[:, 1])
    Y2 = replacevalues(Y, p.N, xyz_n[:, 2])
    Z2 = replacevalues(Z, p.N, xyz_n[:, 3])

    # Store values for continuity in gradients
    return NetworkResults(X2,
        Y2,
        Z2,
        q)
end

function target(r::NetworkResults, p::NetworkOptParams)
    norm(p.X - r.X) + norm(p.Y - r.Y) + norm(p.Z - p.Z)
end

function target(r::NetworkResults, target::Matrix{Float64})
    norm(target .- [r.X r.Y r.Z])
end