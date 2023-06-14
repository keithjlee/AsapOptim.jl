"""
    ktruss(E::Float64, A::Float64, L::Float64)

Get the element truss stiffness matrix in the local coordinate system
"""
function ktruss(E::Float64, A::Float64, L::Float64)
    E * A / L * [1 -1; -1 1]
end

function ktruss(E, A, L)
    E * A / L * [1 -1; -1 1]
end

"""
Adjoint w/r/t element variables E, A, L for the local stiffness matrix of a truss element
"""
function ChainRulesCore.rrule(::typeof(ktruss), E::Float64, A::Float64, L::Float64)
    k = ktruss(E, A, L)

    function ktruss_pullback(k̄)

        # ∇E = dot(k̄, (A / L * [1 -1; -1 1]))
        ∇A = dot(k̄, (E / L * [1 -1; -1 1]))
        ∇L = dot(k̄, (- E * A / L^2 * [1 -1; -1 1]))

        return (NoTangent(), NoTangent(), ∇A, ∇L)
    end

    return k, ktruss_pullback
end

function ChainRulesCore.rrule(::typeof(ktruss), E, A, L)
    k = ktruss(E, A, L)

    function ktruss_pullback(k̄)

        # ∇E = dot(k̄, (A / L * [1 -1; -1 1]))
        ∇A = dot(k̄, (E / L * [1 -1; -1 1]))
        ∇L = dot(k̄, (- E * A / L^2 * [1 -1; -1 1]))

        return (NoTangent(), NoTangent(), ∇A, ∇L)
    end

    return k, ktruss_pullback
end

"""
    getlocalks(E::Vector{Float64}, A::Vector{Float64}, L::Vector{Float64})

Broadcast function for vectors of element properties
"""
function getlocalks(E::Vector{Float64}, A::Vector{Float64}, L::Vector{Float64})
    ktruss.(E, A, L)
end

function getlocalks(E, A, L)
    ktruss.(E, A, L)
end

"""
    kglobal(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, E::Float64, A::Float64, id::Vector{Int64})

Get the element truss stiffness matrix in the global coordinate system
"""
function kglobal(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, E::Float64, A::Float64, id::Vector{Int64})

    i1, i2 = id

    veclocal = [X[i2] - X[i1], Y[i2] - Y[i1], Z[i2] - Z[i1]]
    len = norm(veclocal)

    cx, cy, cz = veclocal ./ len
    r = Rtruss(cx, cy, cz)
    kloc = ktruss(E, A, len)

    r' * kloc * r
end