# """
#     updatemodel(p::TrussOptParams, sol::Optimization.SciMLBase.AbstractOptimizationSolution)

# Generate a new structural model from the results of an optimization
# """
# function updatemodel(p::TrussOptParams, sol::Optimization.SciMLBase.AbstractOptimizationSolution)
    
#     #optim variables
#     u = sol.u

#     #final values
#     X = addvalues(p.X, p.indexer.iX, u[p.indexer.iXg])
#     Y = addvalues(p.Y, p.indexer.iY, u[p.indexer.iYg])
#     Z = addvalues(p.Z, p.indexer.iZ, u[p.indexer.iZg])
#     A = replacevalues(p.A, p.indexer.iA, u[p.indexer.iAg])

#     #new model
#     nodes = Vector{TrussNode}()
#     elements = Vector{TrussElement}()
#     loads = Vector{NodeForce}()

#     #new nodes
#     for (node, x, y, z) in zip(p.model.nodes, X, Y, Z)
#         newnode = TrussNode([x, y, z], node.dof)
#         newnode.id = node.id
#         push!(nodes, newnode)
#     end

#     #new elements
#     for (id, e, a, el) in zip(p.nodeids, p.E, A, p.model.elements)
#         newelement = TrussElement(nodes, id, TrussSection(a, e))
#         newelement.id = el.id
#         push!(elements, newelement)
#     end

#     #new loads
#     for load in p.model.loads
#         newload = NodeForce(nodes[load.node.nodeID], load.value)
#         newload.id = load.id
#         push!(loads, newload)
#     end

#     model = TrussModel(nodes, elements, loads)
#     solve!(model)

#     return model

# end

"""
    updatemodel(p::TrussOptParams, u::Vector{Float64})

Generate a new structural model from the results of an optimization
"""
function updatemodel(p::TrussOptParams, u::Vector{Float64})
    
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

function updatenetwork(p::NetworkOptParams, u::Vector{Float64})
    
    #final values
    X = addvalues(p.X, p.indexer.iX, u[p.indexer.iXg])
    Y = addvalues(p.Y, p.indexer.iY, u[p.indexer.iYg])
    Z = addvalues(p.Z, p.indexer.iZ, u[p.indexer.iZg])
    Q = replacevalues(p.q, p.indexer.iQ, u[p.indexer.iQg])

    #new model
    nodes = Vector{FDMnode}()
    elements = Vector{FDMelement}()
    loads = Vector{FDMload}()

    #new nodes
    for (node, x, y, z) in zip(p.network.nodes, X, Y, Z)
        newnode = FDMnode([x, y, z], node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for (q, el) in zip(Q, p.model.elements)
        newelement = FDMelement(nodes, el.iStart, el.iEnd, q)
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.model.loads
        newload = FDMload(nodes[load.node.nodeID], load.value)
        push!(loads, newload)
    end

    network = Network(nodes, elements, loads)
    solve!(network)

    return network

end

# """
# Store the results of an optimization run and generate a new structural model based on results
# """
# struct OptimResults
#     losstrace::Vector{Float64}
#     valtrace::Vector{Vector{Float64}}
#     solvetime::Float64
#     niter::Int64
#     model::Asap.AbstractModel

#     function OptimResults(p::TrussOptParams, sol::Optimization.SciMLBase.AbstractOptimizationSolution)

#         solvedmodel = updatemodel(p, sol)

#         new(copy(p.losstrace),
#             copy(p.valtrace),
#             sol.solve_time,
#             length(p.valtrace),
#             solvedmodel)

#     end
# end

"""
Results of a truss structural analysis. Basis for all objective functions.
"""
mutable struct TrussResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    A::Vector{Float64}
    L::Vector{Float64}
    K::Vector{Matrix{Float64}}
    R::Vector{Matrix{Float64}}
    U::Vector{Float64}
end

mutable struct NetworkResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    Q::Vector{Float64}
end