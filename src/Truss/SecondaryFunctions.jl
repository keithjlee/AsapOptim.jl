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

            #F̄dF/du
            du[ids[i]] += F̄[i] * (Rs[i] * Eks[i])[2,:]

            # for j = 1:3
            #     (dK[i])[j, :] = F̄[i] * (u[ids[i]])[j] .* (Rs[i])[1, :]
            # end

            # for j = 4:6
            #     (dK[i])[j, :] = F̄[i] * (u[ids[i]])[j] .* (Rs[i])[2, :]
            # end

            #F̄dF/dK
            dK[i] = F̄[i] .* kron((u[ids[i]]), (Rs[i])[2,:]')

            

            # F̄dF/dR
            f1, f2, f3, f4, f5, f6 = F̄[i] .* Eks[i] * u[ids[i]]
            dR[i] = [0 f1 0 f2 0 f3; 0 f4 0 f5 0 f6]
            # dR[i] = [f4 f5 f6 0 0 0; 0 0 0 f4 f5 f6]

        end

        return (NoTangent(), du, dK, dR, NoTangent())
    end

    return F, Faxial_pullback
end