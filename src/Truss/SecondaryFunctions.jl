"""
    Faxial(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

Explicit form for extracting the axial forces of an analyzed truss structure
"""
function Faxial(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    #Displacements w/r/t each element
    uEs = [u[id] for id in p.dofids]

    #End forces
    fE = Rs .* Eks .* uEs

    #axial forces
    getindex.(fE, 2)
end

"""
For a single element, where F is the [2×1] vector of end forces in LCS:

dg/du = dg/dF ⋅ dF/du = F̄ ⋅ RK

dg/dK = dg/dF ⋅ dF/dK = F̄ ⋅ uᵀ ⊗ R

dg/dR = dg/dF ⋅ dF/dR = F̄ ⋅ Ku

And since the scalar axial force is always the *second* component of the end force vector, we take the second row of any adjoint based on R, K 
"""
function ChainRulesCore.rrule(::typeof(Faxial), u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
    F = Faxial(u, Eks, Rs, p)

    function Faxial_pullback(F̄)
        du = zero(u)
        dK = zero.(Eks)
        dR = zero.(Rs)
        ids = p.dofids

        for i in eachindex(ids)

            #F̄ ⋅ dF/du
            du[ids[i]] += F̄[i] * (Rs[i] * Eks[i])[2,:]


            #F̄ ⋅ dF/dK
            dK[i] = F̄[i] .* kron((u[ids[i]]'), (Rs[i])[2,:])
            

            # F̄ ⋅ dF/dR
            f1, f2, f3, f4, f5, f6 = F̄[i] .* Eks[i] * u[ids[i]]
            dR[i] = [0 f1 0 f2 0 f3; 0 f4 0 f5 0 f6]

        end

        return (NoTangent(), du, dK, dR, NoTangent())
    end

    return F, Faxial_pullback
end

"""
    σaxial(F::Vector{Float64}, A::Vector{Float64})

Get the absolute value of the axial stress experienced under force F and area A
"""
function σaxial(F::Vector{Float64}, A::Vector{Float64})
    stress = zero(F)

    for i in eachindex(F)
        if F[i] < 0
            stress[i] = -F[i] / A[i]
        else
            stress[i] = F[i] / A[i]
        end
    end

    return stress
end

"""
dg/dF = dg/dσ ⋅ dσ/dF = σ̄ ⋅ F / |F| / A
dg/dA = dg/dσ ⋅ dσ/dA = σ̄ ⋅ -F / A²
"""
function ChainRulesCore.rrule(::typeof(σaxial), F::Vector{Float64}, A::Vector{Float64})
    stress = σaxial(F, A)

    function σaxial_pullback(σ̄)
        s = sign.(F)
        return (NoTangent(), σ̄ .* s ./ A, -σ̄ .* F ./ A.^2)
    end

    return stress, σaxial_pullback
end

"""
    axialforce(t::TrussResults, p::TrussOptParams)
Axial forces in truss structure
"""
function axialforce(t::TrussResults, p::TrussOptParams)
    Faxial(t.U, t.K, t.R, p)
end
export axialforce

"""
    axialstress(t::TrussResults, p::TrussOptParams)
Axial stresses in truss structure
"""
function axialstress(t::TrussResults, p::TrussOptParams)
    F = Faxial(t.U, t.K, t.R, p)
    σaxial(F, t.A)
end
export axialstress


function Flocal(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    #Displacements w/r/t each element
    uEs = [u[id] for id in p.dofids]

    #End forces
    Rs .* Eks .* uEs
end

function ChainRulesCore.rrule(::typeof(Flocal), u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    F = Flocal(u, Eks, Rs, p)

    function Flocal_pullback(F̄)
        du = zero(u)
        dK = zero.(Eks)
        dR = zero.(Rs)
        ids = p.dofids

        for i in eachindex(ids)
            # F̄ ⋅ dF/dR
            ΔR = kron(Eks[i] * u[ids[i]], I(2)) * F̄[i]
            dR[i] = reshape(ΔR, 2, 6)
            # dR[i] = permutedims(reshape(ΔR, 6, 2), (2,1))

            #F̄ ⋅ dF/du
            du[ids[i]] += (Rs[i] * Eks[i])' * F̄[i]

            #F̄ ⋅ dF/dK
            ΔK = kron(u[ids[i]]', Rs[i])' * F̄[i]
            dK[i] = reshape(ΔK, 6, 6)
            # dK[i] = permutedims(reshape(ΔK, 6, 6), (2,1))
        end

        return (NoTangent(), du, dK, dR, NoTangent())

    end

    return F, Flocal_pullback

end

function Faxial2(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
    fvecs = Flocal(u, Eks, Rs, p)

    getindex.(fvecs, 2)
end