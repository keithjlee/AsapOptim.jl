"""
    getevecs(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::TrussOptParams)

Get the [nₑ × 3] matrix where each row is the [x,y,z] vector of an element
"""
function getevecs(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::AbstractOptParams)
    p.C * [X Y Z]
end

"""
    getlengths(XYZ::Matrix{Float64})

Get the [nₑ × 1] vector of element lengths
"""
function getlengths(XYZ::Matrix{Float64})
    norm.(eachrow(XYZ))
end

"""
    getnormalizedevecs(XYZ::Matrix{Float64}, Ls::Vector{Float64})

Get the unit vector representation of elements (local x axis)
"""
function getnormalizedevecs(XYZ::Matrix{Float64}, Ls::Vector{Float64})
    XYZ ./ Ls
end

"""
    ktruss(E::Float64, A::Float64, L::Float64)

Get the element truss stiffness matrix in the local coordinate system
"""
function ktruss(E::Float64, A::Float64, L::Float64)
    E * A / L * [1 -1; -1 1]
end

function getlocalks(E::Vector{Float64}, A::Vector{Float64}, L::Vector{Float64})
    ktruss.(E, A, L)
end

"""
    Rtruss(Cx::Float64, Cy::Float64, Cz::Float64)

Transformation matrix for truss element
"""
function Rtruss(Cx::Float64, Cy::Float64, Cz::Float64)
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

function Rtruss(Cxyz::SubArray)
    Cx, Cy, Cz = Cxyz
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

function Rtruss(XYZn::Matrix{Float64})
    Rtruss.(eachrow(XYZn))
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
    solveU(K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)

Displacement of free DOFs
"""
function solveU(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    id = p.freeids
    cg(K[id, id], p.P[id])
end

"""
    replacevalues(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

Replace the values of `values[indices]` with the values in `newvalues` in a differentiable way. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function replacevalues(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})
    
    v2 = copy(values)
    v2[indices] .= newvalues

    return v2
end

"""
    addvalues(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

Add the values of `increments` to the current values in `values` at `indices`. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function addvalues(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

    v2 = copy(values)
    v2[indices] .+= increments

    return v2
end