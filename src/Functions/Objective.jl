"""
    solve(values::Vector{Float64}, p::TrussOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solvetruss(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Y = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Z = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    A = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

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