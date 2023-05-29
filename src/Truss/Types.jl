abstract type AbstractOptParams end
abstract type AbstractVariable end
abstract type TrussOptVariable end
abstract type AbstractIndexer end

mutable struct SpatialVariable <: TrussOptVariable
    i::Int64 #index of node, e.g. X[i] is the spatial variable
    val::Float64 #value
    lb::Float64 #lower bound of variable
    ub::Float64 #upper bound of variable
    axis::Symbol #which spatial coordinate?
    iglobal::Int64 # position in the vector of active design variables

    function SpatialVariable(nodeindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, validaxes)

        return new(nodeindex, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::Asap.AbstractNode, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, validaxes)

        return new(node.nodeID, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::Asap.AbstractNode, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, validaxes)

        value = node.position[axis2ind[axis]]

        return new(node.nodeID, value, lowerbound, upperbound, axis)
    end
end

mutable struct AreaVariable <: TrussOptVariable
    i::Int64 #index of element, e.g. A[i] is the area variable
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Int64

    function AreaVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
        return new(elementindex, value, lowerbound, upperbound)
    end

    function AreaVariable(element::Asap.AbstractElement, value::Float64, lowerbound::Float64, upperbound::Float64)
        return new(element.elementID, value, lowerbound, upperbound)
    end

    function AreaVariable(element::Asap.AbstractElement, lowerbound::Float64, upperbound::Float64)
        return new(element.elementID, element.section.A, lowerbound, upperbound)
    end
end

mutable struct CoupledVariable <: AbstractVariable
    i::Int64
    referencevariable::TrussOptVariable

    function CoupledVariable(node::TrussNode, ref::SpatialVariable)
        return new(node.nodeID, ref)
    end

    function CoupledVariable(element::TrussElement, ref::AreaVariable)
        return new(element.elementID, ref)
    end
end

const TrussVariable = Union{SpatialVariable, AreaVariable, CoupledVariable}

"""
    TrussOptIndexer

Translation layer between active design variables and their indices in the global structural property vectors
"""
mutable struct TrussOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
end


function populate!(indexer::TrussOptIndexer, var::SpatialVariable)
    field_local, field_global = axis2field[var.axis]

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), var.iglobal)

end

function populate!(indexer::TrussOptIndexer, var::AreaVariable)
    push!(getfield(indexer, :iA), var.i)
    push!(getfield(indexer, :iAg), var.iglobal)

end

function populate!(indexer::TrussOptIndexer, var::CoupledVariable)
    if typeof(var.referencevariable) == SpatialVariable
        field_local, field_global = axis2field[var.referencevariable.axis]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
    else
        push!(getfield(indexer, :iA), var.i)
        push!(getfield(indexer, :iAg), var.referencevariable.iglobal)
    end
end

"""
    TrussOptIndexer(vars::Vector{TrussVariable})

Generate the index translation layer between model parameters and design variables
"""
function TrussOptIndexer(vars::Vector{TrussVariable})
    indexer = TrussOptIndexer(Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}())

    for var in vars
        populate!(indexer, var)
    end
    
    return indexer
end

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
    losstrace::Vector{Float64} #
    valtrace::Vector{Vector{Float64}}
    # Lstore::Vector{Float64} #current element lengths
    # Rstore::Vector{Matrix{Float64}} #current element transformation matrices
    # Kstore::Vector{Matrix{Float64}} #current element stiffness matrices (GCS)


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
            if typeof(var) <: TrussOptVariable
                var.iglobal  = i
                i += 1
                push!(vals, var.val)
                push!(lowerbounds, var.lb)
                push!(upperbounds, var.ub)
            end
        end

        #generate indexer between design variables and truss parameters
        indexer = TrussOptIndexer(variables)

        #information about the structural model
        nodeids = getproperty.(model.elements, :nodeIDs)
        dofids = getproperty.(model.elements, :globalID)
        P = model.P
        C = Asap.connectivity(model)
        freeids = model.freeDOFs
        inzs = allinz(model)
        cp = model.S.colptr
        rv = model.S.rowval
        nnz = length(model.S.nzval)

        #generate a truss optimization problem
        return new(model, 
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
            model.nDOFs,
            Vector{Float64}(),
            Vector{Vector{Float64}}()
            )

    end
end

"""
    updatemodel(p::TrussOptParams, sol::Optimization.SciMLBase.AbstractOptimizationSolution)

Generate a new structural model from the results of an optimization
"""
function updatemodel(p::TrussOptParams, sol::Optimization.SciMLBase.AbstractOptimizationSolution)
    
    #optim variables
    u = sol.u

    #final values
    X = addvalues(p.X, p.indexer.iX, u[p.indexer.iXg])
    Y = addvalues(p.Y, p.indexer.iY, u[p.indexer.iYg])
    Z = addvalues(p.Z, p.indexer.iZ, u[p.indexer.iZg])
    A = replacevalues(p.A, p.indexer.iA, u[p.indexer.iAg])

    #new model
    nodes = Vector{TrussNode}()
    elements = Vector{TrussElement}()
    loads = Vector{NodeForce}()

    #new nodes
    for (node, x, y, z) in zip(p.model.nodes, X, Y, Z)
        newnode = TrussNode([x, y, z], node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for (id, e, a, el) in zip(p.nodeids, p.E, A, p.model.elements)
        newelement = TrussElement(nodes, id, TrussSection(a, e))
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.model.loads
        newload = NodeForce(nodes[load.node.nodeID], load.value)
        newload.id = load.id
        push!(loads, newload)
    end

    model = TrussModel(nodes, elements, loads)
    solve!(model)

    return model

end

"""
Store the results of an optimization run and generate a new structural model based on results
"""
struct OptimResults
    losstrace::Vector{Float64}
    valtrace::Vector{Vector{Float64}}
    solvetime::Float64
    niter::Int64
    model::Asap.AbstractModel

    function OptimResults(p::TrussOptParams, sol::Optimization.SciMLBase.AbstractOptimizationSolution)

        solvedmodel = updatemodel(p, sol)

        return new(copy(p.losstrace),
            copy(p.valtrace),
            sol.solve_time,
            length(p.valtrace),
            solvedmodel)

    end
end