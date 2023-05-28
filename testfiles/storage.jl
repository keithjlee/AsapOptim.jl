p = TrussOptParams(model, vars);
C = p.C;

#design variables
vals = p.values
indexer = p.indexer

#update values
Xnew = addvalues(p.X, indexer.iX, vals[indexer.iXg])
Ynew = addvalues(p.Y, indexer.iY, vals[indexer.iYg])
Znew = addvalues(p.Z, indexer.iZ, vals[indexer.iZg])
Anew = replacevalues(p.A, indexer.iA, vals[indexer.iAg])
Enew = p.E
#matrix of nodal positions
xyz = [Xnew Ynew Znew]

#functions and adjoints
begin
    function getevecs(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::TrussOptParams)
        p.C * [X Y Z]
    end

    function ChainRulesCore.rrule(::typeof(getevecs), X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::TrussOptParams)
        v = getevecs(X, Y, Z, p)

        function getevecs_pullback(v̄)
            
            # dv = v̄ * ones(3,3)
            
            dv = p.C' * v̄

            return (NoTangent(), dv[:,1], dv[:,2], dv[:,3], NoTangent())
            
        end

        return v, getevecs_pullback
    end

    function getlengths(XYZ::Matrix{Float64})
        norm.(eachrow(XYZ))
    end

    function ChainRulesCore.rrule(::typeof(getlengths), XYZ::Matrix{Float64})
        l = getlengths(XYZ)

        function l_pullback(l̄)
            dl = l̄ ./ l .* XYZ 
            
            return (NoTangent(), dl)
        end

        return l, l_pullback
    end

    function getnormalizedevecs(XYZ::Matrix{Float64}, Ls::Vector{Float64})
        XYZ ./ Ls
    end

    function ChainRulesCore.rrule(::typeof(getnormalizedevecs), XYZ::Matrix{Float64}, Ls::Vector{Float64})
        XYZn = getnormalizedevecs(XYZ, Ls)

        function getnormv_pullback(v̄)
            dxyz = v̄ ./ Ls .* ones(size(XYZn)...)

            dL = sum.(eachrow(-XYZn .* v̄ ./ Ls))

            return (NoTangent(), dxyz, dL)
        end

        return XYZn, getnormv_pullback
    end

    function getglobalks(rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})
        transpose.(rs) .* ks .* rs
    end

    function ChainRulesCore.rrule(::typeof(getglobalks), rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})
        
        kgs = getglobalks(rs, ks)

        function kgs_pullback(k̄)
            dr = (2 .* ks .* rs) .* k̄
            dk = rs .* k̄ .* transpose.(rs)

            return (NoTangent(), dr, dk)
        end

        return kgs, kgs_pullback
    end
end
#element vectors

function objMat(values::Vector{Float64}, p::TrussOptParams)
    indexer = p.indexer

    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])
    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])

    #element vectors
    elementvecs = p.C * [Xnew Ynew Znew]

    #element lengths
    elementlengths = norm.(eachrow(elementvecs))

    #normalized vecs
    elementvecsnormalized = elementvecs ./ elementlengths

    #Γ
    rotmats = AsapOptim.Rtruss.(eachrow(elementvecsnormalized))

    #kloc
    klocs = AsapOptim.klocal.(p.E, Anew, elementlengths)

    #kglob
    kglobs = transpose.(rotmats) .* klocs .* rotmats

    #global stiffness matrix
    K = AsapOptim.assembleglobalK(kglobs, p)

    #solve displacement
    U = solveU(K, p)

    return U' * p.P[p.freeids]
end

@time o = objMat(vals, problem);
@time g = Zygote.gradient(var -> objMat(var, p), vals)[1];

function objMat2(values::Vector{Float64}, p::TrussOptParams)


    indexer = p.indexer
    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])
    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])

    #element vectors
    vₑ = getevecs(Xnew, Ynew, Znew, p)

    #element lengths
    Lₑ = getlengths(vₑ)

    #normalized vecs
    nₑ = getnormalizedevecs(vₑ, Lₑ)

    #transformation matrices
    Γₑ = AsapOptim.Rtruss.(eachrow(nₑ))

    #local stiffness matrices
    kₑ = AsapOptim.klocal.(p.E, Anew, Lₑ)

    #global stiffness matrices
    Kₑ = getglobalks(Γₑ, kₑ)

    #global stiffness matrix
    K = AsapOptim.assembleglobalK(Kₑ, p)

    #solve displacement
    U = solveU(K, p)

    return U' * p.P[p.freeids]
end

func2 = Optimization.OptimizationFunction(objMat2, Optimization.AutoZygote())
prob2 = Optimization.OptimizationProblem(func2, vals, p;
    lb = p.lb,
    ub = p.ub)

function cb2(vals::Vector{Float64}, loss::Float64)
    push!(p.losstrace, loss)
    push!(p.valtrace, deepcopy(vals))
    false
end

cleartrace!(p)
sol2 = Optimization.solve(prob2, NLopt.LD_LBFGS();
    callback = cb2,
    reltol = 1e-3)

res2 = OptimResults(p, sol2)


