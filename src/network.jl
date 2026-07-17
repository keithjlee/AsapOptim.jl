"""
    QVariable <: AbstractVariable

A force-density variable for FDM network optimization: the design entry
REPLACES the element's force density `q` [force/length] (absolute
semantics, like `AreaVariable`).

    QVariable(element::FDMelement, value, lb, ub)
"""
struct QVariable <: AbstractVariable
    element::Asap.FDMelement
    value::Float64
    lb::Float64
    ub::Float64

    function QVariable(element::Asap.FDMelement, value::Real, lb::Real, ub::Real)
        @assert 0 < lb <= value <= ub "need 0 < lb ≤ value ≤ ub (tension form-finding)"
        return new(element, value, lb, ub)
    end
end

"""
    NetworkOptParams

Compiled force-density optimization problem over an `Asap.Network`
(process!-ed): constant connectivity/loads/anchors plus the design-vector
scatter for [`QVariable`](@ref)s (grouping via `CoupledVariable` with an
`FDMelement` target).

The FDM solve is linear in the free coordinates:
`(Cnᵀ diag(q) Cn) · xyz_free = Pn − Cnᵀ diag(q) Cf · xyz_fixed` — evaluated
purely by [`solve_network`](@ref), differentiable end-to-end w.r.t. `q`.
"""
struct NetworkOptParams
    network::Asap.Network
    values::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    q0::Vector{Float64}
    Sq::SparseMatrixCSC{Float64,Int}
    qmask::Vector{Bool}
    Cn::SparseMatrixCSC{Float64,Int}
    Cf::SparseMatrixCSC{Float64,Int}
    Pn::Matrix{Float64}
    xyz_f::Matrix{Float64}
    embed_free::SparseMatrixCSC{Float64,Int}    # n_nodes × n_free row scatter
    embed_fixed::SparseMatrixCSC{Float64,Int}   # n_nodes × n_fixed
    solver::Any                                 # linear-solver backend (see OptParams)
end

"""
    NetworkOptParams(network::Asap.Network, variables::Vector{<:AbstractVariable})

Compile a network and its [`QVariable`](@ref)s (plus `CoupledVariable`
groupings targeting `FDMelement`s) into a [`NetworkOptParams`](@ref):
one design slot per independent variable, a sparse force-density scatter
`Sq`, and constant copies of the connectivity, loads, and anchor positions.
Processes the network first if needed.
"""
function NetworkOptParams(network::Asap.Network, variables::Vector{<:AbstractVariable};
    solver = nothing)
    network.cache === nothing && Asap.process!(network)
    cache = network.cache
    cache.mixed && error("NetworkOptParams does not support per-axis (mixed) node fixity yet — " *
                         "the differentiable path solves all three coordinates against one partition")
    Asap._refresh_state!(cache, network)   # current anchor positions + loads

    slot = Dict{AbstractVariable,Int}()
    values = Float64[]
    lb = Float64[]
    ub = Float64[]
    for v in variables
        v isa CoupledVariable && continue
        v isa QVariable || error("NetworkOptParams takes QVariables (and their couplings); got $(typeof(v))")
        push!(values, v.value)
        push!(lb, v.lb)
        push!(ub, v.ub)
        slot[v] = length(values)
    end

    nel = length(network.elements)
    sI = Int[]
    sJ = Int[]
    sV = Float64[]
    qmask = falses(nel)
    register!(el, j, factor) = begin
        i = el.index
        qmask[i] && error("element $i has two force-density variables")
        qmask[i] = true
        push!(sI, i)
        push!(sJ, j)
        push!(sV, factor)
    end
    for v in variables
        if v isa QVariable
            register!(v.element, slot[v], 1.0)
        elseif v isa CoupledVariable
            v.target isa Asap.FDMelement ||
                error("network couplings must target FDMelements")
            register!(v.target, slot[v.parent], v.factor)
        end
    end

    n = length(network.nodes)
    Nfree = cache.N
    Ffix = cache.F
    ef = sparse(Nfree, 1:length(Nfree), ones(length(Nfree)), n, length(Nfree))
    eF = sparse(Ffix, 1:length(Ffix), ones(length(Ffix)), n, length(Ffix))

    return NetworkOptParams(network, values, lb, ub,
        Float64[el.q for el in network.elements], sparse(sI, sJ, sV, nel, length(values)),
        collect(qmask),
        SparseMatrixCSC{Float64,Int}(cache.C[:, Nfree]), SparseMatrixCSC{Float64,Int}(cache.C[:, Ffix]),
        Matrix{Float64}(cache.P[Nfree, :]), Matrix{Float64}(cache.xyz[Ffix, :]),
        ef, eF, solver)
end

"""
    NetworkOptResults

Outputs of a force-density design evaluation: full nodal coordinates
(`X`, `Y`, `Z` columns of the form-found geometry), the evaluated force
densities `Q`, and member lengths `L`. Member forces are `Q .* L`.
"""
struct NetworkOptResults{TX,TQ,TL}
    X::TX
    Y::TX
    Z::TX
    Q::TQ
    L::TL
end

"""
    solve_network(x, p::NetworkOptParams) -> NetworkOptResults

Evaluate force densities: assemble the FDM system purely and solve for the
free-node coordinates (Asap's `solve_free` multi-RHS adjoint carries the
gradients). Differentiable w.r.t. `x` end-to-end.
"""
function solve_network(x::AbstractVector, p::NetworkOptParams)
    q = p.q0 .* .!p.qmask .+ (p.Sq * x) .* p.qmask

    D = Diagonal(q)
    K = sparse(p.Cn' * D * p.Cn)
    rhs = p.Pn - p.Cn' * (D * (p.Cf * p.xyz_f))

    # default keeps the 2-arg call (Enzyme's imported rules match it)
    xyz_free = p.solver === nothing ? Asap.solve_free(K, rhs) :
               Asap.solve_free(p.solver, K, rhs)
    xyz = p.embed_free * xyz_free + p.embed_fixed * p.xyz_f

    L = _network_lengths(xyz, p)
    return NetworkOptResults(xyz[:, 1], xyz[:, 2], xyz[:, 3], q, L)
end

"""
    _network_lengths(xyz, p::NetworkOptParams) -> Vector

Per-element member lengths of the form-found geometry: row norms of
`C · xyz`, where `C` is the network's signed incidence matrix.
"""
function _network_lengths(xyz, p::NetworkOptParams)
    vx = p.network.cache.C * xyz
    return vec(sqrt.(sum(abs2, vx; dims=2)))
end

"""
    member_forces(res::NetworkOptResults) -> Vector

FDM member forces [force]: `q · L`, tension-positive.
"""
member_forces(res::NetworkOptResults) = res.Q .* res.L
