# method one, direct kglobal
X = problem.X; Y = problem.Y; Z = problem.Z
E = problem.E; A = problem.A
ids = problem.params.nodeids

i = 40



function kgexplicit(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, E::Float64, A::Float64, id::Vector{Int64})

    x1, x2 = X[id]
    y1, y2 = Y[id]
    z1, z2 = Z[id]

    len = L(x1, x2, y1, y2, z1, z2)

    lvec = AsapToolkit.localvector(x1, x2, y1, y2, z1, z2)
    len = L(x1, x2, y1, y2, z1, z2)

    cx, cy, cz = lvec ./ len

    r = Rtruss(cx, cy, cz)
    k = AsapToolkit.klocal(E, A, len)

    r' * k * r
end;

@time kg1 = kglobal(X, Y, Z, E[i], A[i], ids[i]);
@time kg2 = kgexplicit(X, Y, Z, E[i], A[i], ids[i]);

A0 = A[i]

@time g1 = Zygote.gradient(var -> norm(kglobal(X, Y, Z, E[i], var, ids[i])), A0)[1];
@time g2 = Zygote.gradient(var -> norm(kgexplicit(X, Y, Z, E[i], var, ids[i])), A0)[1];

using BenchmarkTools
@btime [kglobal(X, Y, Z, e, a, id) for (e,a,id) in zip(E, A, ids)];
@btime [kgexplicit(X, Y, Z, e, a, id) for (e,a,id) in zip(E, A, ids)];