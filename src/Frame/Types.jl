abstract type FrameOptVariable <: AbstractVariable end

const validproperties = [:A, :Iz, :Iy, :Izz, :Iyy]

const field2field = Dict(:X => (:iX, :iXg),
    :x => (:iX, :iXg),
    :Y => (:iY, :iYg),
    :y => (:iY, :iYg),
    :Z => (:iZ, :iZg),
    :z => (:iZ, :iZg), 
    :A => (:iA, :iAg),
    :Iz => (:iIz, :iIzg),
    :Izz => (:iIz, :iIzg),
    :Iy => (:iIy, :iIyg),
    :Iyy => (:iIy, :iIyg))

mutable struct InternalVariable <: FrameOptVariable
    i::Int64 #index of element, e.g. A[i] is the area variable
    val::Float64
    lb::Float64
    ub::Float64
    property::Symbol
    iglobal::Int64

    function InternalVariable(element::Asap.Element, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol)
        @assert in(property, keys(validproperties))
        new(element.elementID, value, lowerbound, upperbound, property)
    end

    function InternalVariable(element::Asap.Element, lowerbound::Float64, upperbound::Float64, property::Symbol)
        @assert in(property, keys(validproperties))
        new(element.elementID, getfield(element.section, property), lowerbound, upperbound, property)
    end
end

mutable struct FrameOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
    iIz::Vector{Int64}
    iIzg::Vector{Int64}
    iIy::Vector{Int64}
    iIyg::Vector{Int64}
end

function populate!(indexer::FrameOptIndexer, var::SpatialVariable)
    field_local, field_global = field2field[var.axis]

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), var.iglobal)

end

function populate!(indexer::FrameOptIndexer, var::AreaVariable)
    push!(getfield(indexer, :iA), var.i)
    push!(getfield(indexer, :iAg), var.iglobal)
end

function populate!(indexer::FrameOptIndexer, var::InternalVariable)
    field_local, field_global = property2field[var.property]

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), var.iglobal)
end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable)
    if typeof(var.referencevariable) == SpatialVariable
        field_local, field_global = field2field[var.referencevariable.axis]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
    else if typeof(var.referencevariable) == InternalVariable
        field_local, field_global = property2field[var.property]
    
        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
    else
        push!(getfield(indexer, :iA), var.i)
    push!(getfield(indexer, :iAg), var.iglobal)
    end
end

const FrameVariable = Union{SpatialVariable, AreaVariable, InternalVariable, CoupledVariable}
function FrameOptIndexer(vars::Vector{FrameVariable})
    indexer = TrussOptIndexer(
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}()
        )

    for var in vars
        populate!(indexer, var)
    end
    
    indexer
end

mutable struct FrameOptParams <: AbstractOptParams
    model::Model #the reference truss model for optimization
    values::Vector{Float64} #design variables
    indexer::FrameOptIndexer #pointers to design variables and full variables
    variables::Vector{FrameVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    E::Vector{Float64} #all element young's modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    G::Vector{Float64} #shear modulii
    Iz::Vector{Float64} #strong axis moment of inertia
    Iy::Vector{Float64} #weak axis moment of inertia
    J::Vector{Float64} #Torsion constant
    Ψ::Vector{Float64} #LCS orientation
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


    function FrameOptParams(model::Model, variables::Vector{TrussVariable})
        @assert model.processed

        #extract global parameters
        xyz = Asap.nodePositions(model)
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]

        #material parameters
        sections = getproperty.(model.elements, :section)
        E = getproperty.(sections, :E)
        A = getproperty.(sections, :A)
        G = getproperty.(sections, :G)
        Iz = getproperty.(sections, :Izz)
        Iy = getproperty.(sections, :Iyy)
        J = getproperty.(sections, :J)
        Ψ = getproperty.(model.elements, :Ψ)

        #assign global id to variables
        vals = Vector{Float64}()
        lowerbounds = Vector{Float64}()
        upperbounds = Vector{Float64}()

        #assign an index to all unique variables, collect value and bounds
        i = 1
        for var in variables
            if typeof(var) != CoupledVariable
                var.iglobal  = i
                i += 1
                push!(vals, var.val)
                push!(lowerbounds, var.lb)
                push!(upperbounds, var.ub)
            end
        end

        #generate indexer between design variables and truss parameters
        indexer = FrameOptIndexer(variables)

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
            G,
            Iz,
            Iy,
            J,
            Ψ,
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
            Vector{Vector{Float64}}(),
            Vector{Float64}(),
            Vector{Matrix{Float64}}(),
            Vector{Matrix{Float64}}()
            )

    end
end
