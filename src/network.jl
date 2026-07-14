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
end

function NetworkOptParams(network::Asap.Network, variables::Vector{<:AbstractVariable})
    network.processed || Asap.process!(network)

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
        i = el.elementID
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
    Nfree = network.N
    Ffix = network.F
    ef = sparse(Nfree, 1:length(Nfree), ones(length(Nfree)), n, length(Nfree))
    eF = sparse(Ffix, 1:length(Ffix), ones(length(Ffix)), n, length(Ffix))

    return NetworkOptParams(network, values, lb, ub,
        Vector{Float64}(network.q), sparse(sI, sJ, sV, nel, length(values)),
        collect(qmask),
        SparseMatrixCSC{Float64,Int}(network.Cn), SparseMatrixCSC{Float64,Int}(network.Cf),
        Matrix{Float64}(network.Pn), Matrix{Float64}(network.xyz[Ffix, :]),
        ef, eF)
end

"""
    NetworkResults

Outputs of a force-density design evaluation: full nodal coordinates
(`X`, `Y`, `Z` columns of the form-found geometry), the evaluated force
densities `Q`, and member lengths `L`. Member forces are `Q .* L`.
"""
struct NetworkResults{TX,TQ,TL}
    X::TX
    Y::TX
    Z::TX
    Q::TQ
    L::TL
end

"""
    solve_network(x, p::NetworkOptParams) -> NetworkResults

Evaluate force densities: assemble the FDM system purely and solve for the
free-node coordinates (Asap's `solve_free` multi-RHS adjoint carries the
gradients). Differentiable w.r.t. `x` end-to-end.
"""
function solve_network(x::AbstractVector, p::NetworkOptParams)
    qvar = p.Sq * x
    q = [p.qmask[i] ? qvar[i] : p.q0[i] for i in eachindex(p.qmask)]

    D = Diagonal(q)
    K = sparse(p.Cn' * D * p.Cn)
    rhs = p.Pn - p.Cn' * (D * (p.Cf * p.xyz_f))

    xyz_free = Asap.solve_free(K, rhs)
    xyz = p.embed_free * xyz_free + p.embed_fixed * p.xyz_f

    L = _network_lengths(xyz, p)
    return NetworkResults(xyz[:, 1], xyz[:, 2], xyz[:, 3], q, L)
end

function _network_lengths(xyz, p::NetworkOptParams)
    C = p.network.C
    vx = C * xyz
    return [sqrt(vx[i, 1]^2 + vx[i, 2]^2 + vx[i, 3]^2) for i in 1:size(vx, 1)]
end

"""
    member_forces(res::NetworkResults) -> Vector

FDM member forces [force]: `q · L`, tension-positive.
"""
member_forces(res::NetworkResults) = res.Q .* res.L
