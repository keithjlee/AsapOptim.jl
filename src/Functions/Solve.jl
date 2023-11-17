"""
    solve_u(K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)

Displacement of free DOFs
"""
function solve_u(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    id = p.freeids
    cg(K[id, id], p.P[id])
end

function ChainRulesCore.rrule(::typeof(solve_u), K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    u = solve_u(K, p)

    function solve_u_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        dKΔ = cg(K[p.freeids, p.freeids], ū)

        #assign to proper indices
        dudK[p.freeids, p.freeids] .= kron(u', dKΔ)

        return NoTangent(), -dudK, NoTangent()
    end

    return u, solve_u_pullback
end