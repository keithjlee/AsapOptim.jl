"""
    solveU(K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)

Displacement of free DOFs
"""
function solveU(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    id = p.freeids
    cg(K[id, id], p.P[id])
end

function solveU(K, p::AbstractOptParams)
    id = p.freeids
    cg(K[id, id], p.P[id])
end

"""
u = inv(K) * P

if obj = f(u), then the gradient of obj with respect to an independent variable x is achieved through the chain rule:

dObj/dx = df/du ⋅ du/dK ⋅ ... = ū ⋅ dy/dK ⋅ ...

For this rule, we are concerned with finding du/dK, or the [ndof × ndof] matrix of sensitivites that we can propagate backwards to the final objective.

Given df/du = [ndof × 1] = ū is the gradient of the objective function with respect to displacements u, the sensitivity is:

du/dK = - uᵀ ⊗ K⁻¹
df/dK = du/dK ū = - (uᵀ ⊗ K⁻¹)ū

Can be rearranged such that:
ΔK = K⁻¹ū

df/dK = -uᵀ ⊗ ΔK

Which is an [ndof × ndof] matrix where:

Columnᵢ = uᵢ .* ΔK
"""
function ChainRulesCore.rrule(::typeof(solveU), K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    u = solveU(K, p)

    function solveU_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        dKΔ = cg(K[p.freeids, p.freeids], ū)

        #assign to proper indices
        dudK[p.freeids, p.freeids] .= kron(u', dKΔ)

        return NoTangent(), -dudK, NoTangent()
    end

    return u, solveU_pullback
end

function ChainRulesCore.rrule(::typeof(solveU), K, p::AbstractOptParams)
    u = solveU(K, p)

    function solveU_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        dKΔ = cg(K[p.freeids, p.freeids], ū)

        #assign to proper indices
        dudK[p.freeids, p.freeids] .= kron(u', dKΔ)

        return NoTangent(), -dudK, NoTangent()
    end

    return u, solveU_pullback
end