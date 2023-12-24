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
    K::SparseMatrixCSC{Float64, Int64} #[nfreedof × nfreedof] global stiffness matrix
    P::Vector{Float64} # [nfreedof × 1] load vector
    u::Vector{Float64} #[ndof × 1] displacement vector
    # lincache::LinearSolve.LinearCache

    function TrussOptProblem()
    end

    function TrussOptProblem(model::TrussModel)

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

        K = copy(model.S[model.freeDOFs, model.freeDOFs])
        u = copy(model.u)
        P = copy(model.P[model.freeDOFs])

        # lp = LinearProblem(copy(K[model.freeDOFs, model.freeDOFs]), copy(model.P[model.freeDOFs]))
        # lc = init(lp, alg)

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
            P,
            u
            # lc
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
    shadow.P = zero(problem.P)
    shadow.K = explicit_zero(problem.K)
    shadow.u = zero(shadow.u)
    
    # fill!(nonzeros(shadow.lincache.A), 0.)
    # fill!(shadow.lincache.u, 0.)

    return shadow

end

# step by step
function update_values!(x::Vector{Float64}, prob::TrussOptProblem, params::AsapOptim.AbstractOptParams)

    #update values
    params.indexer.activeX && (prob.XYZ[params.indexer.iX, 1] += x[params.indexer.iXg])
    params.indexer.activeY && (prob.XYZ[params.indexer.iY, 2] += x[params.indexer.iYg])
    params.indexer.activeZ && (prob.XYZ[params.indexer.iZ, 3] += x[params.indexer.iZg])
    params.indexer.activeA && (prob.A[params.indexer.iA] = x[params.indexer.iAg])

    nothing

end

function element_update!(prob::TrussOptProblem, params::TrussOptParamsNonalloc)
    #get local stiffness matrix, transformation matrix, and global stiffness matrix
    @simd for i in eachindex(prob.ke)
        prob.Γ[i][1, 1:3] = prob.n[i,:]
        prob.Γ[i][2, 4:6] = prob.n[i,:]
        prob.ke[i] = (params.E[i] * prob.A[i] / prob.L[i] * [1 -1; -1 1])
        prob.Ke[i] = prob.Γ[i]' * prob.ke[i] * prob.Γ[i]
    end
end

function assemble_K!(prob::TrussOptProblem, params::TrussOptParamsNonalloc)
    fill!(nonzeros(prob.K), 0.)
    for i in eachindex(prob.Ke)
        prob.K.nzval[params.inzs[i]] += prob.Ke[i][params.i_k_active[i], params.i_k_active[i]][:]
    end
end

function solve_u!(prob::TrussOptProblem, params::TrussOptParamsNonalloc)
    prob.u[params.freeids] = prob.K \ prob.P
end

function solve_u!(u::Vector{Float64}, K::SparseMatrixCSC{Float64, Int64}, P::Vector{Float64})
    u[params.freeids] = K \ P
end


function solve_u(prob::TrussOptProblem)

    lp = LinearProblem(prob.K, prob.P)
    ls = init(lp, KrylovJL_CG())
    LinearSolve.solve!(ls)

    copy(ls.u)
end

function solve_norm(prob::TrussOptProblem)

    LP = LinearProblem(prob.K, prob.P)
    ls = LinearSolve.solve(LP, LUFactorization())

    norm(ls.u)

end

function solve_norm2(prob::TrussOptProblem, params::TrussOptParamsNonalloc)

    @inbounds lp = LinearProblem(prob.K[params.freeids, params.freeids], params.P[params.freeids])
    ls = LinearSolve.solve(lp, KLUFactorization())

    norm(ls.u)
end

function solve_norm3(K::SparseMatrixCSC{Float64, Int64}, params::TrussOptParamsNonalloc)

    lp = LinearProblem(copy(K), params.P[params.freeids])
    sol = LinearSolve.solve!(lp, KLUFactorization())

    norm(sol.u)
end

function solve_norm4(K::SparseMatrixCSC{Float64, Int64}, P::Vector{Float64})

    lp = LinearProblem(K, P)
    # ls = init(lp, KLUFactorization())
    sol = LinearSolve.solve(lp, LUFactorization())

    norm(sol.u)
end

export nonalloc
function nonalloc(x::Vector{Float64}, prob::AbstractTrussOptProblem, params::TrussOptParamsNonalloc)

    #update values
    update_values!(x, prob, params)

    #element vectors
    prob.v = params.C * prob.XYZ

    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        @views prob.L[i] = norm(prob.v[i, :])
        @views prob.n[i, :] = prob.v[i, :] / prob.L[i]
    end

    #update element-wise stiffness/transformation information
    element_update!(prob, params)

    # sum(norm.(prob.Ke))

    #assemble global K
    assemble_K!(prob, params)

    norm(prob.K)

    #solve and norm
    # solve_norm4(prob.K, params.P)
    # solve_norm(prob)

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

    # norm(K[p.freeids, p.freeids])
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
