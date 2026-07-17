"""
    OptParams

The compiled optimization problem: everything constant during optimization,
prepared once so each design evaluation is a handful of sparse products and
one pure Asap solve.

Built from a model and its design variables:

    OptParams(model, variables)

(`TrussOptParams` and `FrameOptParams` are accepted aliases — the v1.0 core
is unified, so one parameter type serves both.)

# Fields (all constant per optimization)
- `model::Model`: the processed reference model (its cache holds the frozen
  sparsity pattern, DOF partition, and load vectors)
- `values`, `lb`, `ub::Vector{Float64}`: the design vector's start and bounds
- `X0::Matrix{Float64}`: base node positions (3 × n_nodes)
- `Sx::SparseMatrixCSC`: position scatter — `X(x) = X0 + reshape(Sx·x, 3, :)`
  (spatial variables are ADDITIVE perturbations)
- `A0::Vector{Float64}`: per-element base areas (elements without an area
  variable keep these)
- `Sa::SparseMatrixCSC`, `amask::Vector{Bool}`: area scatter — for elements
  with a variable, `A(x) = Sa·x` REPLACES the area (absolute semantics)
- `F::Vector{Float64}`: the free-DOF load vector (loads are constant data
  in the differentiable path)
- `i1`, `i2::Vector{Int}`: per-element node indices; `base_sections`: the
  reference sections — plain-data mirrors so the differentiable path never
  reads mutable structs (a Zygote tangent-accumulation catastrophe at scale)

# Design-evaluation contract
`x` is the vector an optimizer manipulates; [`solve_structure`](@ref)`(x, p)`
gives results, gradients flow through Zygote (or any ChainRules-aware
engine) with no additional rules.
"""
struct OptParams
    model::Model{Float64}
    values::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    X0::Matrix{Float64}
    Sx::SparseMatrixCSC{Float64,Int}
    A0::Vector{Float64}
    Sa::SparseMatrixCSC{Float64,Int}
    amask::Vector{Bool}
    F::Vector{Float64}
    i1::Vector{Int}
    i2::Vector{Int}
    base_sections::Vector{Any}
    Evec::Vector{Float64}            # per-element Young's moduli (constant)
    # joint-stiffness variables: design slot (0 = none) and factor per
    # element end — rotational springs ky = kz replaced by factor·x[slot]
    jslot1::Vector{Int}
    jfac1::Vector{Float64}
    jslot2::Vector{Int}
    jfac2::Vector{Float64}
    base_ends::Vector{Any}           # per-element reference EndConditions
end

const TrussOptParams = OptParams
const FrameOptParams = OptParams

"""
    _axis_component(axis::Symbol) -> Int

Row offset (1/2/3) of a global axis (`:X`/`:Y`/`:Z`) within a node's column
of the 3 × n position matrix.
"""
_axis_component(axis::Symbol) = axis === :X ? 1 : axis === :Y ? 2 : 3

"""
    OptParams(model::Model, variables::Vector{<:AbstractVariable})

Compile a processed (or processable) `Asap.Model` and its design variables
into an [`OptParams`](@ref): assign one design-vector slot per independent
variable, build the sparse position/area scatter maps (couplings become
extra scattered entries with their factors), record the joint-variable
slots, and snapshot all constant data (base positions, sections, loads) as
plain arrays.

Errors early on ill-posed declarations: two variables on one element's
area or joint, couplings whose parent is not among the independents,
target/parent type mismatches, and joint variables on elements that carry
element loads (their fixed-end forces would not track the design).
"""

function OptParams(model::Model{Float64}, variables::Vector{<:AbstractVariable})
    model.cache === nothing && process!(model)
    cache = model.cache
    assemble_loads!(cache, model)

    # assign one design-vector slot per independent variable
    slot = Dict{AbstractVariable,Int}()
    values = Float64[]
    lb = Float64[]
    ub = Float64[]
    for v in variables
        v isa CoupledVariable && continue
        push!(values, v.value)
        push!(lb, v.lb)
        push!(ub, v.ub)
        slot[v] = length(values)
    end
    nx = length(values)

    nel = length(model.elements)
    nnodes = length(model.nodes)

    # scatter maps
    sxI = Int[]
    sxJ = Int[]
    sxV = Float64[]
    saI = Int[]
    saJ = Int[]
    saV = Float64[]
    amask = falses(nel)

    register_spatial!(node, j, factor, axis) = begin
        push!(sxI, 3 * (node.index - 1) + _axis_component(axis))
        push!(sxJ, j)
        push!(sxV, factor)
    end
    register_area!(el, j, factor) = begin
        amask[el.index] && error("element $(el.index) (:$(el.id)) has two area variables")
        amask[el.index] = true
        push!(saI, el.index)
        push!(saJ, j)
        push!(saV, factor)
    end

    jslot1 = zeros(Int, nel)
    jfac1 = ones(nel)
    jslot2 = zeros(Int, nel)
    jfac2 = ones(nel)
    register_joint!(el, position, j, factor) = begin
        el isa FrameElement || error("joint variables require a FrameElement")
        if position in (:start, :both)
            jslot1[el.index] == 0 || error("element $(el.index) start joint has two variables")
            jslot1[el.index] = j
            jfac1[el.index] = factor
        end
        if position in (:end, :both)
            jslot2[el.index] == 0 || error("element $(el.index) end joint has two variables")
            jslot2[el.index] = j
            jfac2[el.index] = factor
        end
    end

    for v in variables
        if v isa SpatialVariable
            register_spatial!(v.node, slot[v], 1.0, v.axis)
        elseif v isa AreaVariable
            register_area!(v.element, slot[v], 1.0)
        elseif v isa JointVariable
            register_joint!(v.element, v.position, slot[v], 1.0)
        elseif v isa CoupledVariable
            haskey(slot, v.parent) ||
                error("coupled variable's parent is not among the independent variables")
            j = slot[v.parent]
            if v.parent isa SpatialVariable
                v.target isa Node ||
                    error("a variable coupled to a SpatialVariable must target a Node")
                register_spatial!(v.target, j, v.factor, v.parent.axis)
            elseif v.parent isa JointVariable
                tgt, pos = v.target isa Tuple ? v.target : (v.target, v.parent.position)
                tgt isa FrameElement ||
                    error("a variable coupled to a JointVariable must target a FrameElement (or (element, position) tuple)")
                register_joint!(tgt, pos, j, v.factor)
            else
                v.target isa Union{FrameElement,TrussElement} ||
                    error("a variable coupled to an AreaVariable must target an element")
                register_area!(v.target, j, v.factor)
            end
        end
    end

    X0 = zeros(3, nnodes)
    for (i, n) in enumerate(model.nodes)
        X0[:, i] = n.position
    end
    A0 = [el.section isa Section ? el.section.A : 0.0 for el in model.elements]

    F = cache.P .- cache.Pf
    i1 = [el.nodeStart.index for el in model.elements]
    i2 = [el.nodeEnd.index for el in model.elements]
    # element-load fixed-end forces depend on end conditions, but the pure
    # path treats loads as constant data — a joint variable on an element
    # that carries element loads would silently use stale FEFs. Guard it.
    for load in model.loads
        if load isa ElementLoad && load.element isa FrameElement
            i = load.element.index
            (jslot1[i] != 0 || jslot2[i] != 0) && error(
                "element $(i) (:$(load.element.id)) carries an element load AND a joint " *
                "variable — its fixed-end forces would not track the design. Apply the " *
                "load as NodeForces, or split the member at the load point.")
        end
    end

    base_sections = Any[el.section for el in model.elements]
    Evec = [s isa Section ? s.material.E : 0.0 for s in base_sections]
    base_ends = Any[el isa Union{FrameElement,VariableElement} ? el.ends :
                    EndConditions(:fixedfixed) for el in model.elements]

    return OptParams(model, values, lb, ub, X0,
        sparse(sxI, sxJ, sxV, 3 * nnodes, nx),
        A0, sparse(saI, saJ, saV, nel, nx), collect(amask), F,
        i1, i2, base_sections, Evec,
        jslot1, jfac1, jslot2, jfac2, base_ends)
end
