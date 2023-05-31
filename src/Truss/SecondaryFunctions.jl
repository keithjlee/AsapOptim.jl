function Faxial(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    #Displacements w/r/t each element
    uEs = [u[id] for id in p.dofids]

    #End forces
    fE = Rs .* Eks .* uEs

    #axial forces
    getindex.(fE, 2)
end

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

function ChainRulesCore.rrule(::typeof(σaxial), F::Vector{Float64}, A::Vector{Float64})
    stress = σaxial(F, A)

    function σaxial_pullback(σ̄)
        s = sign.(F)
        return (NoTangent(), σ̄ .* s ./ A, -σ̄ .* F ./ A.^2)
    end

    return stress, σaxial_pullback
end