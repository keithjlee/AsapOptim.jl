abstract type AbstractTrussOptProblem end

mutable struct TrussOptProblem <: AbstractTrussOptProblem
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
    # lincache::LinearSolve.LinearCache

    function TrussOptProblem()
    end

    function TrussOptProblem(model::TrussModel; alg = KLUFactorization())

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

        # lp = LinearProblem(K[model.freeDOFs, model.freeDOFs], model.P[model.freeDOFs])
        # ls = init(lp, alg)

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
            # ls
        )


    end

end

export shadow
function shadow(problem::TrussOptProblem)

    shadow = deepcopy(problem)

    shadow.XYZ = zero(problem.XYZ)
    shadow.A = zero(problem.A)
    shadow.v = zero(problem.v)
    shadow.L = zero(problem.L)
    shadow.n = zero(problem.n)
    shadow.Γ = zero.(problem.Γ)
    shadow.ke = zero.(problem.ke)
    shadow.Ke = zero.(problem.Ke)
    shadow.K = explicit_zero(problem.K)
    shadow.u = zero(shadow.u)

    return shadow

end

function element_update!(prob::TrussOptProblem, params::TrussOptParams)
    #get local stiffness matrix, transformation matrix, and global stiffness matrix
    for i in eachindex(prob.ke)
        prob.Γ[i][1, 1:3] = prob.n[i,:]
        prob.Γ[i][2, 4:6] = prob.n[i,:]
        prob.ke[i] .= (params.E[i] * prob.A[i] / prob.L[i] * [1 -1; -1 1])
        prob.Ke[i] .= prob.Γ[i]' * prob.ke[i] * prob.Γ[i]
    end
end

function assemble_K!(prob::TrussOptProblem, params::TrussOptParams)
    fill!(nonzeros(prob.K), 0.)
    for i in eachindex(prob.Ke)
        prob.K.nzval[params.inzs[i]] += prob.Ke[i][:]
        # prob.K[params.dofids[i], params.dofids[i]] += prob.Ke[i]
    end
end

function solve_u!(prob::TrussOptProblem, params::TrussOptParams)
    prob.u[params.freeids] = prob.K[params.freeids, params.freeids] \ params.P[params.freeids]
end

function solve_u!(u::Vector{Float64}, K::SparseMatrixCSC{Float64, Int64}, params::TrussOptParams)
    u[params.freeids] = K[params.freeids, params.freeids] \ params.P[params.freeids]
end


function solve_u(prob::TrussOptProblem, params::TrussOptParams)

    lp = LinearProblem(prob.K[params.freeids, params.freeids], params.P[params.freeids])
    ls = init(lp, KLUFactorization())
    LinearSolve.solve!(ls)

    copy(ls.u)
end

export nonalloc
function nonalloc(x::Vector{Float64}, prob::AbstractTrussOptProblem, params::TrussOptParams)

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

    #update element-wise stiffness/transformation information
    element_update!(prob, params)

    # sum(norm.(prob.Ke))

    #assemble global K
    assemble_K!(prob, params)

    # norm(prob.K)
    # # solve_u!(prob, params)
    # # solve_u!(prob.u, prob.K, params)

    lp = LinearProblem(prob.K[params.freeids, params.freeids], params.P[params.freeids])
    ls = LinearSolve.solve(lp, LUFactorization())
    norm(ls.u)

    # prob.u[params.freeids] = ls.u

    # # prob.u[params.freeids] .= prob.K[params.freeids, params.freeids] \ params.P[params.freeids]
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

    # sum(norm.(Kₑ))
    # # K
    K = assemble_global_K(Kₑ, p)

    # norm(K)
    # # norm(K)
    # # K⁻¹P
    u = solve_u_direct(K, p)
    norm(u)

    # # U
    # U = replace_values(zeros(p.n), p.freeids, u)
    # norm(U)

    # norm(U)
    # U' * p.P
end
