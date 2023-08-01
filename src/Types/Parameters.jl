"""
    TrussOptParams(model::TrussModel, variables::Vector{TrussVariable})

Contains all information and fields necessary for optimization.
"""
mutable struct TrussOptParams <: AbstractOptParams
    model::TrussModel #the reference truss model for optimization
    values::Vector{Float64} #design variables
    indexer::TrussOptIndexer #pointers to design variables and full variables
    variables::Vector{TrussVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    E::Vector{Float64} #all element young's modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    P::Vector{Float64} # External load vector
    C::SparseMatrixCSC{Int64, Int64} #connectivity matrix
    lb::Vector{Float64} #lower bounds of variables
    ub::Vector{Float64} #upper bounds of variables
    cp::Vector{Int64} #S.colptr
    rv::Vector{Int64} #S.rowval
    nnz::Int64 #length(S.nzval)
    inzs::Vector{Vector{Int64}} # Indices of elemental K in global S.nzval
    freeids::Vector{Int64} # [DofFree1, DofFree2,...]
    nodeids::Vector{Vector{Int64}} # [[iNodeStart, iNodeEnd] for element in elements]
    dofids::Vector{Vector{Int64}} # [[dofStartNode..., dofEndNode...] for element in elements]
    n::Int64 #total number of DOFs

    function TrussOptParams(model::TrussModel, variables::Vector{TrussVariable})
        @assert model.processed

        #extract global parameters
        xyz = Asap.nodePositions(model)
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]
        E = getproperty.(getproperty.(model.elements, :section), :E)
        A = getproperty.(getproperty.(model.elements, :section), :A)

        #assign global id to variables
        vals = Vector{Float64}()
        lowerbounds = Vector{Float64}()
        upperbounds = Vector{Float64}()

        #assign an index to all unique variables, collect value and bounds
        i = 1
        for var in variables
            if typeof(var) <: Union{SpatialVariable, AreaVariable}
                var.iglobal  = i
                i += 1
                push!(vals, var.val)
                push!(lowerbounds, var.lb)
                push!(upperbounds, var.ub)
            end
        end

        #generate indexer between design variables and truss parameters
        indexer = TrussOptIndexer(variables)

        #topology
        nodeids = getproperty.(model.elements, :nodeIDs)
        dofids = getproperty.(model.elements, :globalID)
        C = Asap.connectivity(model)
        freeids = model.freeDOFs

        #external load
        P = model.P

        #sparsity pattern of K
        inzs = allinz(model)
        cp = model.S.colptr
        rv = model.S.rowval
        nnz = length(model.S.nzval)

        #generate a truss optimization problem
        new(model, 
            vals, 
            indexer, 
            variables, 
            X, 
            Y, 
            Z, 
            E, 
            A, 
            P,
            C,
            lowerbounds, 
            upperbounds,
            cp,
            rv,
            nnz,
            inzs,
            freeids,
            nodeids,
            dofids,
            model.nDOFs
            )

    end
end

"""
    NetworkOptParams(model::Network, variables::Vector{NetworkVariable})

Contains all information and fields necessary for optimization.
"""
mutable struct NetworkOptParams <: AbstractOptParams
    network::Network #the reference truss model for optimization
    values::Vector{Float64} #design variables
    indexer::NetworkOptIndexer #pointers to design variables and full variables
    variables::Vector{NetworkVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    q::Vector{Float64} #all element young's modulii |n_element|
    Pn::Matrix{Float64}
    C::SparseMatrixCSC{Int64, Int64} #connectivity matrix
    Cn::SparseMatrixCSC{Int64, Int64}
    Cf::SparseMatrixCSC{Int64, Int64}
    lb::Vector{Float64} #lower bounds of variables
    ub::Vector{Float64} #upper bounds of variables
    N::Vector{Int64}
    F::Vector{Int64}

    function NetworkOptParams(network::Network, variables::Vector{NetworkVariable})
        @assert network.processed "network must be processed"

        #extract global parameters
        xyz = network.xyz
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]
        q = network.q

        #assign global id to variables
        vals = Vector{Float64}()
        lowerbounds = Vector{Float64}()
        upperbounds = Vector{Float64}()

        #assign an index to all unique variables, collect value and bounds
        i = 1
        for var in variables
            if typeof(var) <: Union{SpatialVariable, QVariable}
                var.iglobal  = i
                i += 1
                push!(vals, var.val)
                push!(lowerbounds, var.lb)
                push!(upperbounds, var.ub)
            end
        end

        #generate indexer between design variables and truss parameters
        indexer = NetworkOptIndexer(variables)

        #generate a truss optimization problem
        new(network,
            vals,
            indexer,
            variables,
            X,
            Y,
            Z,
            q,
            network.Pn,
            network.C,
            network.Cn,
            network.Cf,
            lowerbounds,
            upperbounds,
            network.N,
            network.F
            )

    end
end