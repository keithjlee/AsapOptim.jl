"""
    OptResults

Outputs of one design evaluation ([`solve_structure`](@ref)), carrying
everything downstream objectives and constraints consume — all plain data,
so gradients of any scalar of these fields flow back to the design vector.

# Fields
- `U::Vector`: full-space displacements (6 slots per node — index a node's
  vertical displacement as `U[6(i−1)+3]`)
- `X::Matrix`: evaluated node positions (3 × n_nodes)
- `A::Vector`: evaluated per-element areas
- `L::Vector`: evaluated element lengths
- `sections::AbstractVector`: the evaluated per-element sections (feed to
  further pure computations if needed)
- `EA::AbstractVector`: evaluated axial rigidities (plain vector — cheap
  under reverse AD)
"""
struct OptResults{TU,TX,TA,TL,TS,TE}
    U::TU
    X::TX
    A::TA
    L::TL
    sections::TS
    EA::TE
end

# design vector → (positions, areas, sections, state); the shared pure core
function _design_state(x::AbstractVector, p::OptParams)
    X = p.X0 + reshape(p.Sx * x, 3, :)
    Avar = p.Sa * x
    A = [p.amask[i] ? Avar[i] : p.A0[i] for i in eachindex(p.amask)]
    sections = map(eachindex(p.amask)) do i
        s = p.base_sections[i]::Section{Float64}
        p.amask[i] ? Section(s.material, A[i], s.Ix, s.Iy, s.J) : s
    end
    EAvec = p.Evec .* A
    return X, A, sections, EAvec
end

_element_lengths(X, p::OptParams) =
    map(eachindex(p.i1)) do k
        i, j = p.i1[k], p.i2[k]
        sqrt((X[1, j] - X[1, i])^2 + (X[2, j] - X[2, i])^2 + (X[3, j] - X[3, i])^2)
    end

"""
    solve_structure(x, p::OptParams) -> OptResults

Evaluate a design vector: scatter positions and areas, run Asap's pure
solve on the resulting `ModelState`, and package the results. Fully
differentiable — `Zygote.gradient(x -> f(solve_structure(x, p)), x0)`
needs nothing else.

`solve_truss` and `solve_frame` are aliases (the unified core needs no
distinction).
"""
function solve_structure(x::AbstractVector, p::OptParams)
    X, A, sections, EAvec = _design_state(x, p)
    state = ModelState{Float64}(X, sections, EAvec)
    U = Asap.solve(p.model, state)
    return OptResults(U, X, A, _element_lengths(X, p), sections, EAvec)
end

const solve_truss = solve_structure
const solve_frame = solve_structure

"""
    axial_force(res::OptResults, p::OptParams) -> Vector

Per-element axial force [force], tension-positive: `EA/L` times the axial
elongation extracted from the global displacements. Pure and
differentiable; exact for truss elements and for the axial action of frame
elements with rigid axial connections.
"""
function axial_force(res::OptResults, p::OptParams)
    map(eachindex(p.i1)) do i
        ni, nj = p.i1[i], p.i2[i]
        L = res.L[i]
        nx = (res.X[1, nj] - res.X[1, ni]) / L
        ny = (res.X[2, nj] - res.X[2, ni]) / L
        nz = (res.X[3, nj] - res.X[3, ni]) / L
        du1 = res.U[6*(nj-1)+1] - res.U[6*(ni-1)+1]
        du2 = res.U[6*(nj-1)+2] - res.U[6*(ni-1)+2]
        du3 = res.U[6*(nj-1)+3] - res.U[6*(ni-1)+3]
        res.EA[i] / L * (nx * du1 + ny * du2 + nz * du3)
    end
end

"""
    axial_stress(res::OptResults, p::OptParams) -> Vector

Per-element axial stress [force/length²]: axial force over evaluated area.
"""
axial_stress(res::OptResults, p::OptParams) = axial_force(res, p) ./ res.A

"""
    compliance(res::OptResults, p::OptParams) -> Real

External work `Fᵀu` of the evaluated design — the canonical smooth
stiffness objective.
"""
compliance(res::OptResults, p::OptParams) = dot(p.F, res.U)

"""
    GeometricProperties(x, p::OptParams)

Geometry-only evaluation of a design — positions, lengths, areas — with NO
structural solve. Use for objectives like material volume
(`dot(geo.A, geo.L)`) where paying for a solve would be waste.

# Fields
- `X::Matrix`, `L::Vector`, `A::Vector`
"""
struct GeometricProperties{TX,TL,TA}
    X::TX
    L::TL
    A::TA
end

function GeometricProperties(x::AbstractVector, p::OptParams)
    # geometry-only: skip section-struct construction entirely (its
    # constructor pullbacks would dominate this otherwise trivial path)
    X = p.X0 + reshape(p.Sx * x, 3, :)
    Avar = p.Sa * x
    A = [p.amask[i] ? Avar[i] : p.A0[i] for i in eachindex(p.amask)]
    return GeometricProperties(X, _element_lengths(X, p), A)
end
