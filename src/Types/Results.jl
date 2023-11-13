"""
    updatemodel(p::TrussOptParams, u::Vector{Float64})

Generate a new structural model from the results of an optimization
"""
function updatemodel(p::TrussOptParams, u::Vector{Float64})
    
    #final values
    X = addvalues(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = addvalues(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = addvalues(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)
    A = replacevalues(p.A, p.indexer.iA, u[p.indexer.iAg] .* p.indexer.fA)

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
    updatenetwork(p::TrussOptParams, u::Vector{Float64})

Generate a new FDM network from the results of an optimization
"""
function updatenetwork(p::NetworkOptParams, u::Vector{Float64})
    
    #final values
    X = addvalues(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = addvalues(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = addvalues(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)
    Q = replacevalues(p.q, p.indexer.iQ, u[p.indexer.iQg] .* p.indexer.fQ)

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
    for (q, el) in zip(Q, p.network.elements)
        newelement = FDMelement(nodes, el.iStart, el.iEnd, q)
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.network.loads
        newload = FDMload(nodes[load.point.nodeID], load.force)
        push!(loads, newload)
    end

    network = Network(nodes, elements, loads)
    solve!(network)

    return network

end

"""
    TrussResults

Results of a differentiable structural analysis of a truss model.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- A: Area of all elements
- L: Length of all elements
- K: Elemental stiffness matrices in GCS
- R: Elemental transformation matrices
- U: Displacement vector of all nodes
"""
struct TrussResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    A::Vector{Float64}
    L::Vector{Float64}
    K::Vector{Matrix{Float64}}
    R::Vector{Matrix{Float64}}
    U::Vector{Float64}
end

"""
    NetworkResults

Results of a differentiable FDM analysis of a network.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- Q: force densities of all elements
"""
struct NetworkResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    Q::Vector{Float64}
end

"""
    GeometricProperties

Get the primary and secondary geometric properties of a truss structure in a differentiable format.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- evecs: [nel × 3] matrix where row_i = nodes[iend].position - nodes[istart].position
- A: Area of all elements
- L: Length of all elements
- 
"""
struct GeometricProperties
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    evecs::Matrix{Float64}
    A::Vector{Float64}
    L::Vector{Float64}

    function GeometricProperties(design_variables::Vector{Float64}, opt_params::TrussOptParams)

        #populate values
        X = addvalues(opt_params.X, opt_params.indexer.iX, design_variables[opt_params.indexer.iXg] .* opt_params.indexer.fX)
        Y = addvalues(opt_params.Y, opt_params.indexer.iY, design_variables[opt_params.indexer.iYg] .* opt_params.indexer.fY)
        Z = addvalues(opt_params.Z, opt_params.indexer.iZ, design_variables[opt_params.indexer.iZg] .* opt_params.indexer.fZ)
        A = replacevalues(opt_params.A, opt_params.indexer.iA, design_variables[opt_params.indexer.iAg] .* opt_params.indexer.fA)

        # [nₑₗ × 3] matrix where row i is the vector representation of element i, from the start node to the end node; ||vecₑ|| = Lₑ
        evecs = getevecs(X, Y, Z, p)
        L = getlengths(evecs)

        return new(
            X,
            Y,
            Z,
            evecs,
            A,
            L
        )

    end
end
