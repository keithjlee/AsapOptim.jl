struct FrameOptParams <: AbstractOptParams
    model::Model #the reference  model for optimization
    values::Vector{Float64} #design variables
    indexer::FrameOptIndexer #pointers to design variables and full variables
    variables::Vector{AbstractVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    Ψ::Vector{Float64} #all roll angles |n_element|
    E::Vector{Float64} #all element young's modulii |n_element|
    G::Vector{Float64} #all element shear modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    Ix::Vector{Float64} #all element X-X moment of inertia (nominal strong) |n_element|
    Iy::Vector{Float64} #all element Y-Y moment of inertia (nominal weak) |n_element|
    J::Vector{Float64} #all element torsion constant |n_element|
    P::Vector{Float64} # External load vector
    Pf::Vector{Float64} # Fixed end force vector
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

    function FrameOptParams(model::Model, variables::Vector{T}) where T <: FrameVariable

        #assert model is processed
        model.processed || (Asap.process!(model))

        #global parameters
        
        #spatial positions
        xyz = node_positions(model)
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]
        Ψ = getproperty.(model.elements, :Ψ)

        #Material properties
        sections = getproperty.(model.elements, :section)

        E = getproperty.(sections, :E)
        G = getproperty.(sections, :G)

        #Geometric properties
        A = getproperty.(sections, :A)
        Ix = getproperty.(sections, :Ix)
        Iy = getproperty.(sections, :Iy)
        J = getproperty.(sections, :J)

        #Variable processing
        vals, lowerbounds, upperbounds = process_variables!(variables)

        #Indexer
        indexer = FrameOptIndexer(variables)

        #topology
        nodeids = getproperty.(model.elements, :nodeIDs)
        dofids = getproperty.(model.elements, :globalID)
        C = Asap.connectivity(model)
        freeids = model.freeDOFs

        #loads
        P = model.P
        Pf = model.Pf

        #sparsity pattern of K
        inzs = all_inz(model)
        cp = model.S.colptr
        rv = model.S.rowval
        nnz = length(model.S.nzval)


        new(
            model,
            vals,
            indexer,
            variables,
            X,
            Y,
            Z,
            Ψ,
            E,
            G,
            A,
            Ix,
            Iy,
            J,
            P,
            Pf,
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

export FrameOptParams2
struct FrameOptParams2 <: AbstractOptParams
    model::Model #the reference  model for optimization
    # values::Vector{Float64} #design variables
    # indexer::FrameOptIndexer #pointers to design variables and full variables
    # variables::Vector{AbstractVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    Ψ::Vector{Float64} #all roll angles |n_element|
    E::Vector{Float64} #all element young's modulii |n_element|
    G::Vector{Float64} #all element shear modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    Ix::Vector{Float64} #all element X-X moment of inertia (nominal strong) |n_element|
    Iy::Vector{Float64} #all element Y-Y moment of inertia (nominal weak) |n_element|
    J::Vector{Float64} #all element torsion constant |n_element|
    # P::Vector{Float64} # External load vector
    # Pf::Vector{Float64} # Fixed end force vector
    # C::SparseMatrixCSC{Int64, Int64} #connectivity matrix
    # lb::Vector{Float64} #lower bounds of variables
    # ub::Vector{Float64} #upper bounds of variables
    # cp::Vector{Int64} #S.colptr
    # rv::Vector{Int64} #S.rowval
    # nnz::Int64 #length(S.nzval)
    # inzs::Vector{Vector{Int64}} # Indices of elemental K in global S.nzval
    # freeids::Vector{Int64} # [DofFree1, DofFree2,...]
    # nodeids::Vector{Vector{Int64}} # [[iNodeStart, iNodeEnd] for element in elements]
    # dofids::Vector{Vector{Int64}} # [[dofStartNode..., dofEndNode...] for element in elements]
    # n::Int64 #total number of DOFs

    function FrameOptParams2(model::Model, variables::Vector{T}) where T <: FrameVariable

        #assert model is processed
        model.processed || (Asap.process!(model))

        #global parameters
        
        #spatial positions
        xyz = node_positions(model)
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]
        Ψ = getproperty.(model.elements, :Ψ)

        #Material properties
        sections = getproperty.(model.elements, :section)

        E = getproperty.(sections, :E)
        G = getproperty.(sections, :G)

        #Geometric properties
        A = getproperty.(sections, :A)
        Ix = getproperty.(sections, :Ix)
        Iy = getproperty.(sections, :Iy)
        J = getproperty.(sections, :J)

        #Variable processing
        # vals, lowerbounds, upperbounds = process_variables!(variables)

        # #Indexer
        # indexer = FrameOptIndexer(variables)

        # #topology
        # nodeids = getproperty.(model.elements, :nodeIDs)
        # dofids = getproperty.(model.elements, :globalID)
        # C = Asap.connectivity(model)
        # freeids = model.freeDOFs

        # #loads
        # P = model.P
        # Pf = model.Pf

        # #sparsity pattern of K
        # inzs = all_inz(model)
        # cp = model.S.colptr
        # rv = model.S.rowval
        # nnz = length(model.S.nzval)


        new(
            model,
            # vals,
            # indexer,
            # variables,
            X,
            Y,
            Z,
            Ψ,
            E,
            G,
            A,
            Ix,
            Iy,
            J,
            # P,
            # Pf,
            # C,
            # lowerbounds,
            # upperbounds,
            # cp,
            # rv,
            # nnz,
            # inzs,
            # freeids,
            # nodeids,
            # dofids,
            # model.nDOFs
        )


    end

end