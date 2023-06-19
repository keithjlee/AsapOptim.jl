const dRdx_truss = [1. 0. 0. 0. 0. 0.; 0. 0. 0. 1. 0. 0.]
const dRdy_truss = [0. 1. 0. 0. 0. 0.; 0. 0. 0. 0. 1. 0.]
const dRdz_truss = [0. 0. 1. 0. 0. 0.; 0. 0. 0. 0. 0. 1.]

"""
    Rtruss(Cx::Float64, Cy::Float64, Cz::Float64)

Transformation matrix for truss element
"""
function Rtruss(Cx::Float64, Cy::Float64, Cz::Float64)
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

"""
Adjoint for global transformation matrix
"""
function ChainRulesCore.rrule(::typeof(Rtruss), Cx::Float64, Cy::Float64, Cz::Float64)
    R = Rtruss(Cx, Cy, Cz)

    function Rtruss_pullback(R̄)
        return (NoTangent(),
            dot(R̄, dRdx_truss),
            dot(R̄, dRdy_truss),
            dot(R̄, dRdz_truss))
    end

    return R, Rtruss_pullback
end

"""
    Rtruss(Cxyz::SubArray)

Transformation of sliced view of matrix of element vectors
"""
function Rtruss(Cxyz::SubArray)
    Cx, Cy, Cz = Cxyz
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

function ChainRulesCore.rrule(::typeof(Rtruss), Cxyz::SubArray)
    R = Rtruss(Cxyz)

    function Rtruss_pullback(R̄)
        return (NoTangent(),
            [dot(R̄, dRdx_truss),
            dot(R̄, dRdy_truss),
            dot(R̄, dRdz_truss)])
    end

    return R, Rtruss_pullback
end

"""
    Rtruss(XYZn::Matrix{Float64})

Get all transformation matrices from a [nₑ × 3] matrix of all normalized element local x vectors
"""
function Rtruss(XYZn::Matrix{Float64})
    Rtruss.(eachrow(XYZn))
end

function ChainRulesCore.rrule(::typeof(Rtruss), XYZn::Matrix{Float64})
    rs = Rtruss(XYZn)

    function Rtruss_pullback(R̄)
        dRdC = zero(XYZn)

        for i in axes(R̄, 1)
            dRdC[i, 1] = dot(R̄[i], dRdx_truss)
            dRdC[i, 2] = dot(R̄[i], dRdy_truss)
            dRdC[i, 3] = dot(R̄[i], dRdz_truss)
        end

        return (NoTangent(), dRdC)
    end

    return rs, Rtruss_pullback
end