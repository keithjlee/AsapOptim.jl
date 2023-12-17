# step by step
function update_values!(x::Vector{Float64}, prob::TrussOptProblem, params::TrussOptParams)

    #update values
    params.indexer.activeX && (prob.XYZ[params.indexer.iX, 1] += x[params.indexer.iXg])
    params.indexer.activeY && (prob.XYZ[params.indexer.iY, 2] += x[params.indexer.iYg])
    params.indexer.activeZ && (prob.XYZ[params.indexer.iZ, 3] += x[params.indexer.iZg])
    params.indexer.activeA && (prob.A[params.indexer.iA] = x[params.indexer.iAg])

    nothing

end

function update_element_vectors!(prob::TrussOptProblem, params::TrussOptParams)
    #element vectors
    prob.v = params.C * prob.XYZ

    nothing
end

function update_element_properties!(prob::TrussOptProblem)

    #element lengths, normalized vectors
    @simd for i in axes(prob.v, 1)
        prob.L[i] = norm(prob.v[i, :])
        prob.n[i, :] = prob.v[i, :] ./ prob.L[i]
    end

    nothing
end

function update_element_matrices!(prob::TrussOptProblem, params::TrussOptParams)
    #get local stiffness matrix, transformation matrix, and global stiffness matrix
    @simd for i in eachindex(prob.ke)
        prob.Γ[i][1, 1:3] = prob.n[i,:]
        prob.Γ[i][2, 4:6] = prob.n[i,:]
        prob.ke[i] = (params.E[i] * prob.A[i] / prob.L[i] * [1 -1; -1 1])
        prob.Ke[i] = prob.Γ[i]' * prob.ke[i] * prob.Γ[i]
    end

    nothing
end

function solvenorm(K::SparseMatrixCSC{Float64, Int64}, P::Vector{Float64}; alg = KLUFactorization())
    LP = LinearProblem(K, P)
    # LS = init(LP, alg)
    LS = LinearSolve.solve(LP, alg)

    norm(LS.u)
end

x0 = params.values ./ 5 .* (1 .+ rand())

using JET
@report_opt update_values!(x0, prob, params)
@report_opt update_element_vectors!(prob, params)
@report_opt update_element_properties!(prob)
@report_opt update_element_matrices!(prob, params)
@report_opt AsapOptim.assemble_K!(prob, params)

K = copy(prob.K[params.freeids, params.freeids])
P = copy(params.P[params.freeids])

@report_opt solvenorm(K, P)

dK = AsapOptim.explicit_zero(K)
dP = zero(P)

function solvenorm2(K::SparseMatrixCSC{Float64, Int64}, P::Vector{Float64})
    LP = LinearProblem(copy(K), P)
    # LS = init(LP, alg)
    LS = LinearSolve.solve(LP, LUFactorization())

    uvals = LS.u
    norm(uvals)
end

function solvenorm_call(K::SparseMatrixCSC{Float64, Int64}, P::Vector{Float64}, ids::Vector{Int64})

    solvenorm2(K[ids, ids], P[ids])
    
end


begin
    K = copy(prob.K[params.freeids, params.freeids])
    P = copy(params.P[params.freeids])

    dK = AsapOptim.explicit_zero(K)
    dP = zero(P)

    Kfull = copy(prob.K)
    dKfull = AsapOptim.explicit_zero(Kfull)
    Pfull = copy(params.P)
    dPfull = zero(Pfull)
    inds = copy(params.freeids)
end

@time solvenorm2(K, P)
@time norm(K\P)
@time solvenorm3(K, lc)
@time solvenorm_call(Kfull, Pfull, inds)

@time autodiff(
    Enzyme.Reverse,
    solvenorm2,
    Duplicated(K, dK),
    Duplicated(P, dP)
)