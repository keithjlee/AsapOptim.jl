abstract type AbstractTrussOptProblem end
mutable struct TrussOptProblem <: AbstractTrussOptProblem
    XYZ::Matrix{Float64} #[nₙ × 3] nodal positions
    A::Vector{Float64} #[nₑ × 1] element areas
    v::Matrix{Float64} #[nₑ × 3] matrix of element vectors
    L::Vector{Float64} #[nₑ × 1] element lengths
    n::Matrix{Float64} #[nₑ × 3] normalized element vectors
    Γ::Array{Float64, 3} #[2 × 6 × nₑ] of transformation matrices
    ke::Array{Float64, 3} #[2 × 2 × nₑ] of stiffness matrices in LCS
    Ke::Array{Float64, 3} #[6 × 6 × nₑ] of stiffness matrices in GCS
    K::SparseMatrixCSC{Float64, Int64} #[ndof × ndof] global stiffness matrix
    u::Vector{Float64} #[ndof × 1] displacement vector

    function TrussOptProblem()
    end

    function TrussOptProblem(model::TrussModel)

        XYZ = node_positions(model)

        #initialize
        A = zeros(model.nElements)
        v = zeros(model.nElements, 3)
        L = zeros(model.nElements)
        n = zeros(model.nElements, 3)

        Γ = Array{Float64, 3}(undef, 2, 6, model.nElements)
        ke = Array{Float64, 3}(undef, 2, 2, model.nElements)
        Ke = Array{Float64, 3}(undef, 6, 6, model.nElements)

        for i in eachindex(model.elements)

            element = model.elements[i]

            A[i] = element.section.A
            L[i] = element.length
            n[i, :] = first(element.LCS)
            v[i, :] = element.length * first(element.LCS)

            Γ[:, :, i] = element.R
            Ke[:, :, i] = element.K
            ke[:, :, i] = element.R * element.K * element.R'

        end

        K = model.S
        u = model.u

        new(
            XYZ,
            A,
            v,
            L,
            n,
            Γ,
            ke,
            Ke,
            K,
            u
        )


    end

end

export TrussOptProblem2
mutable struct TrussOptProblem2 <: AbstractTrussOptProblem
    XYZ::Matrix{Float64} #[nₙ × 3] nodal positions
    A::Vector{Float64} #[nₑ × 1] element areas
    v::Matrix{Float64} #[nₑ × 3] matrix of element vectors
    L::Vector{Float64} #[nₑ × 1] element lengths
    n::Matrix{Float64} #[nₑ × 3] normalized element vectors
    Γ::Vector{MMatrix{2,6,Float64,12}} #[2 × 6 × nₑ] of transformation matrices
    ke::Vector{MMatrix{2,2,Float64,4}} #[2 × 2 × nₑ] of stiffness matrices in LCS
    Ke::Vector{MMatrix{6,6,Float64,36}} #[6 × 6 × nₑ] of stiffness matrices in GCS
    K::SparseMatrixCSC{Float64, Int64} #[ndof × ndof] global stiffness matrix
    u::Vector{Float64} #[ndof × 1] displacement vector

    function TrussOptProblem2()
    end

    function TrussOptProblem2(model::TrussModel)

        XYZ = node_positions(model)

        #initialize
        A = zeros(model.nElements)
        v = zeros(model.nElements, 3)
        L = zeros(model.nElements)
        n = zeros(model.nElements, 3)

        Γ = Vector{MMatrix{2,6,Float64,12}}()
        ke = Vector{MMatrix{2,2,Float64,4}}()
        Ke = Vector{MMatrix{6,6,Float64,36}}()

        for i in eachindex(model.elements)

            element = model.elements[i]

            A[i] = element.section.A
            L[i] = element.length
            n[i, :] = first(element.LCS)
            v[i, :] = element.length * first(element.LCS)

            push!(Γ, MMatrix{2,6}(element.R))
            push!(Ke, MMatrix{6,6}(element.K))
            push!(ke, MMatrix{2,2}(element.R * element.K * element.R'))

        end

        K = model.S
        u = model.u

        new(
            XYZ,
            A,
            v,
            L,
            n,
            Γ,
            ke,
            Ke,
            K,
            u
        )


    end

end

export shadow
function shadow(problem::AbstractTrussOptProblem)

    shadow = deepcopy(problem)

    for field in fieldnames(typeof(problem))

        if typeof(getproperty(shadow, field)) == SparseMatrixCSC{Float64, Int64}
            setproperty!(shadow, field, explicit_zero(getproperty(shadow, field)))
        else
            setproperty!(shadow, field, zero(getproperty(shadow, field)))
        end
    end

    return shadow

end

function solve_u!(prob::AbstractTrussOptProblem, params::TrussOptParams)

    cg!(prob.u[params.freeids], prob.K[params.freeids, params.freeids], params.P[params.freeids])

end

function assemble_K!(K::SparseMatrixCSC{Float64, Int64}, Ke::Array{Float64, 3}, params::TrussOptParams)

    for i in eachindex(K.nzval)
        K.nzval[i] = 0.
    end

    for i in axes(Ke, 3)
        K.nzval[params.inzs[i]] += Ke[:, :, i][:]
        # K[params.dofids[i], params.dofids[i]] += Ke[:, :, i]
    end

    nothing
end

function assemble_K!(prob::TrussOptProblem, params::TrussOptParams)


    fill!(nonzeros(prob.K), 0.)
    @simd for i in axes(prob.Ke, 3)
        prob.K.nzval[params.inzs[i]] += prob.Ke[:, :, i][:]
    end
end

function assemble_K!(prob::TrussOptProblem2, params::TrussOptParams)


    fill!(nonzeros(prob.K), 0.)
    @simd for i in eachindex(prob.Ke)
        @view(prob.K.nzval[params.inzs[i]]) .+= @view(prob.Ke[i][:])
    end

end

export nonalloc
function nonalloc(x::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParams)

    #update values
    params.indexer.activeX && (prob.XYZ[params.indexer.iX, 1] += x[params.indexer.iXg])
    params.indexer.activeY && (prob.XYZ[params.indexer.iY, 2] += x[params.indexer.iYg])
    params.indexer.activeZ && (prob.XYZ[params.indexer.iZ, 3] += x[params.indexer.iZg])
    params.indexer.activeA && (prob.A[params.indexer.iA] = x[params.indexer.iAg])

    #element vectors
    prob.v = params.C * prob.XYZ

    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        prob.L[i] = norm(prob.v[i, :])
        prob.n[i, :] = prob.v[i, :] / prob.L[i]
    end

    #get local stiffness matrix, transformation matrix, and global stiffness matrix
    @simd for i in axes(prob.ke, 3)
        prob.Γ[:, :, i] = r_truss(prob.n[i, :])
        prob.ke[:, :, i] = k_truss(params.E[i], prob.A[i], prob.L[i])
        prob.Ke[:, :, i] = prob.Γ[:, :, i]' * prob.ke[:, :, i] * prob.Γ[:, :, i]

    end

    #assemble global K
    assemble_K!(prob, params)

    norm(prob.K)

    # prob.u[params.freeids] .= prob.K[params.freeids, params.freeids] \ params.P[params.freeids]

    # norm(prob.u)
end

function nonalloc(x::Vector{Float64}, prob::TrussOptProblem2, params::TrussOptParams)

    #update values
    params.indexer.activeX && (prob.XYZ[params.indexer.iX, 1] += x[params.indexer.iXg])
    params.indexer.activeY && (prob.XYZ[params.indexer.iY, 2] += x[params.indexer.iYg])
    params.indexer.activeZ && (prob.XYZ[params.indexer.iZ, 3] += x[params.indexer.iZg])
    params.indexer.activeA && (prob.A[params.indexer.iA] = x[params.indexer.iAg])

    #element vectors
    prob.v = params.C * prob.XYZ

    #element lengths, normalized vectors
    @inbounds @simd for i in axes(prob.v, 1)
        prob.L[i] = norm(prob.v[i, :])
        prob.n[i, :] = prob.v[i, :] / prob.L[i]
    end

    #get local stiffness matrix, transformation matrix, and global stiffness matrix
    @inbounds @simd for i in eachindex(prob.ke)
        prob.Γ[i] .= r_truss_static(prob.n[i, :])
        prob.ke[i] .= (params.E[i] * prob.A[i] / prob.L[i] * [1 -1; -1 1])
        prob.Ke[i] .= prob.Γ[i]' * prob.ke[i] * prob.Γ[i]

    end

    #assemble global K
    assemble_K!(prob, params)

    norm(prob.K)

    # prob.u[params.freeids] .= prob.K[params.freeids, params.freeids] \ params.P[params.freeids]

    # norm(prob.u)
end


export alloc
function alloc(x::Vector{Float64}, p::TrussOptParams)

    #update values
    #populate values
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, x[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, x[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, x[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    A = p.indexer.activeA ? replace_values(p.A, p.indexer.iA, x[p.indexer.iAg] .* p.indexer.fA) : p.A

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

    norm(K)
    # # K⁻¹P
    # u = solve_u(K, p)

    # # U
    # U = replace_values(zeros(p.n), p.freeids, u)

    # norm(U)
    # U' * p.P
end
