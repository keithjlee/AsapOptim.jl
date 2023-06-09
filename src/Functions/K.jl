"""
    getglobalks(rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})

Get a vector of elemental stiffness matrices in GCS given a vector of transformation matrices and a vector of elemental stiffness matrices in LCS
"""
function getglobalks(rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})
    transpose.(rs) .* ks .* rs
end

"""
    assembleglobalK(elementalKs::Vector{Matrix{Float64}}, p::TrussOptParams)

Assemble the global stiffness matrix from a vector of elemental stiffness matrices
"""
function assembleglobalK(elementalKs::Vector{Matrix{Float64}}, p::AbstractOptParams)

    nz = zeros(p.nnz)

    for (k, i) in zip(elementalKs, p.inzs)
        nz[i] .+= vec(k)
    end

    SparseMatrixCSC(p.n, p.n, p.cp, p.rv, nz)
end

"""
The sensitivity of K w/r/t an elemental Kₑ is the proportional stiffness added to K from Kₑ

Output is a vector: [nElements × [nDOFe × nDOFe]] of elemental stiffness matrix sensitivites
"""
function ChainRulesCore.rrule(::typeof(assembleglobalK), Eks::Vector{Matrix{Float64}}, p::AbstractOptParams)
    K = assembleglobalK(Eks, p)

    function K_pullback(K̄)
        dEks = [Matrix(K̄[id, id]) for id in p.dofids]

        return NoTangent(), dEks, NoTangent()
    end

    return K, K_pullback
end