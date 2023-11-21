"""
    solve_truss(values::Vector{Float64}, p::TrussOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_truss(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ)
    A = replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA)

    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    l = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, l)

    # Γ
    Γ = r_truss(n)

    # kₑ
    kₑ = k_truss.(p.E, A, l)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u(K, p)

    # U
    U = replace_values(zeros(p.n), p.freeids, u)

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
    solve_truss_direct(values::Vector{Float64}, p::TrussOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_truss_direct(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ)
    A = replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA)

    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    l = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, l)

    # Γ
    Γ = r_truss(n)

    # kₑ
    kₑ = k_truss.(p.E, A, l)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u_direct(K, p)

    # U
    U = replace_values(zeros(p.n), p.freeids, u)

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
    max_penalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
Penalize values above a threshold
"""
function max_penalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
    Δ = vals .- threshold
    factor * sum(Δ .+ abs.(Δ))
end

"""
    min_penalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
Penalize values below a threshold
"""
function min_penalty(vals::Vector{Float64}, threshold::Float64; factor = 1.)
    Δ = threshold .- vals
    factor * sum(Δ .+ abs.(Δ))
end

"""
    solve_network(values::Vector{Float64}, p::NetworkOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_network(values::Vector{Float64}, p::NetworkOptParams)
    
    #populate values
    X = add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ)
    q = replace_values(p.q, p.indexer.iQ, values[p.indexer.iQg] .* p.indexer.fQ)

    # fixed nodal positions
    xyz_f = [X[p.F] Y[p.F] Z[p.F]]

    # diagonal q matrix
    Q = diagm(q)

    #solve for free positions
    xyz_n = (p.Cn' * Q * p.Cn) \ (p.Pn - p.Cn' * Q * p.Cf * xyz_f)

    X2 = replace_values(X, p.N, xyz_n[:, 1])
    Y2 = replace_values(Y, p.N, xyz_n[:, 2])
    Z2 = replace_values(Z, p.N, xyz_n[:, 3])

    # Store values for continuity in gradients
    return NetworkResults(X2,
        Y2,
        Z2,
        q)
end

function target(r::NetworkResults, p::NetworkOptParams)
    norm([(p.X - r.X) (p.Y - r.Y) (p.Z - p.Z)])
end

function target(r::NetworkResults, target::Matrix{Float64})
    norm(target .- [r.X r.Y r.Z])
end