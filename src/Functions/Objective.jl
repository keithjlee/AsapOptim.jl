"""
    solve_truss(values::Vector{Float64}, p::TrussOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_truss(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    A = p.indexer.activeA ? replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A


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
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    A = p.indexer.activeA ? replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A

    # vₑ 
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
    solve_network(values::Vector{Float64}, p::NetworkOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_network(values::Vector{Float64}, p::NetworkOptParams)
    
    #populate values
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    q = p.indexer.activeQ ? replace_values(p.q, p.indexer.iQ, values[p.indexer.iQg] .* p.indexer.fQ) : p.q

    # fixed nodal positions
    xyz_f = [X[p.F] Y[p.F] Z[p.F]]

    # diagonal q matrix
    Q = sparse(diagm(q))

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

function solve_frame(values::Vector{Float64}, p::FrameOptParams)

    #populate values
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z

    A = p.indexer.activeA ? replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A
    Ix = p.indexer.activeIx ? replace_values(p.Ix, p.indexer.iIx, values[p.indexer.iIxg] .* p.indexer.fIx) : p.Ix
    Iy = p.indexer.activeIy ? replace_values(p.Iy, p.indexer.iIy, values[p.indexer.iIyg] .* p.indexer.fIy) : p.Iy
    J = p.indexer.activeJ ? replace_values(p.J, p.indexer.iJ, values[p.indexer.iJg] .* p.indexer.fJ) : p.J

    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    L = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, L)

    # Γ
    Γ = r_frame(n, p.Ψ)

    # kₑ
    kₑ = k_frame.(p.E, p.G, A, L, Ix, Iy, J)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u(K, p)

    # U
    U = replace_values(zeros(p.n), p.freeids, u)

    return FrameResults(
        X,
        Y,
        Z,
        A,
        Ix,
        Iy,
        J,
        L,
        Kₑ,
        Γ,
        U
    )
end

function solve_frame_direct(values::Vector{Float64}, p::FrameOptParams)

    #populate values
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z

    A = p.indexer.activeA ? replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A
    Ix = p.indexer.activeIx ? replace_values(p.Ix, p.indexer.iIx, values[p.indexer.iIxg] .* p.indexer.fIx) : p.Ix
    Iy = p.indexer.activeIy ? replace_values(p.Iy, p.indexer.iIy, values[p.indexer.iIyg] .* p.indexer.fIy) : p.Iy
    J = p.indexer.activeJ ? replace_values(p.J, p.indexer.iJ, values[p.indexer.iJg] .* p.indexer.fJ) : p.J

    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    L = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, L)

    # Γ
    Γ = r_frame(n, p.Ψ)

    # kₑ
    kₑ = k_frame.(p.E, p.G, A, L, Ix, Iy, J)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u_direct(K, p)

    # U
    U = replace_values(zeros(p.n), p.freeids, u)

    return FrameResults(
        X,
        Y,
        Z,
        A,
        Ix,
        Iy,
        J,
        L,
        Kₑ,
        Γ,
        U
    )
end