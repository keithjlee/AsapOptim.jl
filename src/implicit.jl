#=
Implicit differentiation of the design → solution map.

Differentiating K(x)·u(x) = F (loads design-constant) gives

    du/dxⱼ = −K⁻¹ (∂K/∂xⱼ · u)

so the WHOLE solution tangent matrix costs one assembly, one
factorization, and ONE multi-RHS back-substitution — where chunked
forward-mode AD re-assembles K with Dual numbers once per chunk. The
pseudo-loads (∂K/∂xⱼ)u are element-local: each design variable perturbs
only the kernel inputs of the element(s) it touches.

Generality contract (do not break it): NOTHING in the hot path
enumerates variable types. Stage 1 computes per-element kernel-input
sensitivities — a closed-form fast path where trivial (truss area:
K linear in A), and otherwise a GENERIC one-partial ForwardDiff seed
through the same `truss_stiffness`/`frame_stiffness` kernels the solve
uses, so any variable that enters through kernel arguments (future
section-property variables included) is covered with zero new
derivative code. Stage 2 chains through the compiled scatter maps
(`Sx`, `Sa`, joint slots), which ARE ∂(kernel inputs)/∂x — couplings,
factors, and arbitrary mixing come along structurally.

This is a plain-Float64 evaluation path for optimizer callbacks — it is
not itself differentiable (differentiate `solve_structure` instead).
=#

using ForwardDiff: Dual, Partials, value, partials

struct _ImplicitTag end

_dual1(v::Real) = Dual{_ImplicitTag}(Float64(v), 1.0)

# seed component c of a position with a unit partial
_seed(v::SVector{3,Float64}, c::Int) =
    SVector{3,Dual{_ImplicitTag,Float64,1}}(ntuple(k -> Dual{_ImplicitTag}(v[k], k == c ? 1.0 : 0.0), 3))

"""
    DesignTangents

Solution of [`solution_tangents`](@ref): the design evaluation `res`
(an [`OptResults`](@ref)) plus the full solution tangent matrix
`dU` (n_dofs × n_vars, `dU[:, j] = ∂U/∂xⱼ`). Feed to
[`axial_force_jacobian`](@ref)/[`axial_stress_jacobian`](@ref), or slice
`dU` directly for displacement-constraint rows (e.g. vertical
displacements: `t.dU[3:6:end, :]`).
"""
struct DesignTangents{TR}
    res::TR
    dU::Matrix{Float64}
end

"""
    solution_tangents(x, p::OptParams) -> DesignTangents

Evaluate the design AND the exact solution tangents `∂U/∂x` by implicit
differentiation: one assembly, one factorization (through `p.solver`),
one multi-RHS back-substitution — instead of one Dual re-assembly per
forward-mode chunk. Equal to the ForwardDiff/Zygote Jacobian of the
displacement field to machine precision, at a fraction of the cost
(the per-iteration workhorse for constraint Jacobians in optimization
loops).

Plain `Float64` only — this path is for optimizer callbacks and is not
itself AD-transparent.
"""
function solution_tangents(x::AbstractVector{Float64}, p::OptParams)
    X, A, sections, EAvec, ends = _design_state(x, p)
    cache = p.model.cache
    state = ModelState{Float64}(X, sections, EAvec, ends)
    K = Asap.assemble_K(cache, state)
    fact = Asap._factorize(p.solver, K)

    free = cache.partition.free
    uf = fact \ p.F[free]
    U = cache.free_embed * uf
    res = OptResults(U, X, A, _element_lengths(X, p), sections, EAvec)

    n = length(x)
    nnodes = size(X, 2)
    nel = length(p.i1)

    # stage 2 scatter maps, inverted for stage 1's element loop:
    # node -> [(coord, slot, factor)] and element -> [(slot, factor)]
    node_inputs = [Tuple{Int,Int,Float64}[] for _ in 1:nnodes]
    let rows = rowvals(p.Sx), vals = nonzeros(p.Sx)
        for j in 1:n, ptr in nzrange(p.Sx, j)
            r = rows[ptr]
            push!(node_inputs[(r-1)÷3+1], ((r - 1) % 3 + 1, j, vals[ptr]))
        end
    end
    area_inputs = [Tuple{Int,Float64}[] for _ in 1:nel]
    let rows = rowvals(p.Sa), vals = nonzeros(p.Sa)
        for j in 1:n, ptr in nzrange(p.Sa, j)
            push!(area_inputs[rows[ptr]], (j, vals[ptr]))
        end
    end

    # stage 1: element-local pseudo-loads W[:, j] = (∂K/∂xⱼ)·U
    W = zeros(length(U), n)
    for (e, el) in enumerate(p.model.elements)
        i1, i2 = p.i1[e], p.i2[e]
        touched = !isempty(area_inputs[e]) || !isempty(node_inputs[i1]) ||
                  !isempty(node_inputs[i2]) || p.jslot1[e] != 0 || p.jslot2[e] != 0
        touched || continue
        el isa Union{TrussElement,FrameElement} || error(
            "implicit differentiation does not support design variables touching " *
            "$(typeof(el)) yet — use ForwardDiff through solve_structure instead")

        x1 = SVector{3}(X[1, i1], X[2, i1], X[3, i1])
        x2 = SVector{3}(X[1, i2], X[2, i2], X[3, i2])
        sec = sections[e]::Section{Float64}

        if el isa TrussElement
            dofs = SVector{6,Int}(6i1 - 5, 6i1 - 4, 6i1 - 3, 6i2 - 5, 6i2 - 4, 6i2 - 3)
            ue = SVector{6}(U[d] for d in dofs)

            # area fast path: the truss kernel is LINEAR in A
            if !isempty(area_inputs[e])
                w = (Asap.truss_stiffness(sec, x1, x2) * ue) ./ A[e]
                for (j, fac) in area_inputs[e]
                    @views W[dofs, j] .+= fac .* w
                end
            end
            # spatial: generic one-partial seed through the kernel
            for (xi, ni) in ((1, i1), (2, i2)), (c, j, fac) in node_inputs[ni]
                Kd = xi == 1 ? Asap.truss_stiffness(sec, _seed(x1, c), x2) :
                     Asap.truss_stiffness(sec, x1, _seed(x2, c))
                w = partials.(Kd, 1) * ue
                @views W[dofs, j] .+= fac .* w
            end
        else # FrameElement
            dofs = SVector{12,Int}(ntuple(k -> k <= 6 ? 6(i1 - 1) + k : 6(i2 - 1) + k - 6, 12))
            ue = SVector{12}(U[d] for d in dofs)
            ec = ends === nothing ? el.ends : ends[e]::EndConditions{Float64}
            roll = el.rollangle

            # frame area: seed the section (correct even when semi-rigid
            # axial springs make K nonlinear in A)
            if !isempty(area_inputs[e])
                secd = Section(sec.material, _dual1(A[e]), sec.Ix, sec.Iy, sec.J)
                w = partials.(Asap.frame_stiffness(secd, ec, x1, x2, roll), 1) * ue
                for (j, fac) in area_inputs[e]
                    @views W[dofs, j] .+= fac .* w
                end
            end
            # spatial
            for (xi, ni) in ((1, i1), (2, i2)), (c, j, fac) in node_inputs[ni]
                Kd = xi == 1 ? Asap.frame_stiffness(sec, ec, _seed(x1, c), x2, roll) :
                     Asap.frame_stiffness(sec, ec, x1, _seed(x2, c), roll)
                w = partials.(Kd, 1) * ue
                @views W[dofs, j] .+= fac .* w
            end
            # joint springs: ky = kz share the slot — seed both together
            if p.jslot1[e] != 0
                e1 = ec.e1
                ecd = EndConditions(EndSprings(e1.kx, e1.kt, _dual1(e1.ky), _dual1(e1.kz)), ec.e2)
                w = partials.(Asap.frame_stiffness(sec, ecd, x1, x2, roll), 1) * ue
                @views W[dofs, p.jslot1[e]] .+= p.jfac1[e] .* w
            end
            if p.jslot2[e] != 0
                e2 = ec.e2
                ecd = EndConditions(ec.e1, EndSprings(e2.kx, e2.kt, _dual1(e2.ky), _dual1(e2.kz)))
                w = partials.(Asap.frame_stiffness(sec, ecd, x1, x2, roll), 1) * ue
                @views W[dofs, p.jslot2[e]] .+= p.jfac2[e] .* w
            end
        end
    end

    # one multi-RHS back-substitution for every tangent at once
    dUf = -(fact \ W[free, :])
    return DesignTangents(res, Matrix(cache.free_embed * dUf))
end

"""
    axial_force_jacobian(t::DesignTangents, p::OptParams) -> Matrix

Exact Jacobian `∂N/∂x` (n_elements × n_vars) of the per-element axial
forces, chaining the solution tangents with the closed-form geometry
terms (`N = EA·(ΔX·ΔU)/L²` differentiated in `U`, `X`, and `A`).
"""
function axial_force_jacobian(t::DesignTangents, p::OptParams)
    res = t.res
    X, U, A, L, EA = res.X, res.U, res.A, res.L, res.EA
    nel = length(p.i1)
    ndof = length(U)
    nnodes = size(X, 2)

    ΔX = X * transpose(p.Cinc)                          # 3 × nel
    ΔU = reshape(U, 6, :)[1:3, :] * transpose(p.Cinc)
    s = vec(sum(ΔX .* ΔU; dims=1))
    N = EA .* s ./ (L .^ 2)

    # sparse row operators: row e picks ± the element's 3-vector at its ends
    tI = Int[]; tJ = Int[]; tVx = Float64[]; tVu = Float64[]
    xI = Int[]; xJ = Int[]                              # X-side (3 dofs/node)
    for e in 1:nel
        i1, i2 = p.i1[e], p.i2[e]
        for c in 1:3
            push!(tI, e, e); push!(tJ, 6(i2 - 1) + c, 6(i1 - 1) + c)
            push!(xI, e, e); push!(xJ, 3(i2 - 1) + c, 3(i1 - 1) + c)
            push!(tVx, ΔX[c, e], -ΔX[c, e])
            push!(tVu, ΔU[c, e], -ΔU[c, e])
        end
    end
    S_X6 = sparse(tI, tJ, tVx, nel, ndof)               # ΔX picked at U dofs
    S_Xx = sparse(xI, xJ, tVx, nel, 3nnodes)            # ΔX picked at X dofs
    S_Ux = sparse(xI, xJ, tVu, nel, 3nnodes)            # ΔU picked at X dofs

    cU = EA ./ (L .^ 2)
    J = Diagonal(cU) * (S_X6 * t.dU)                    # U̇ term (dense)
    J .+= Diagonal(cU) * (S_Ux * p.Sx)                  # direct ΔẊ·ΔU term
    dL = Diagonal(1.0 ./ L) * (S_Xx * p.Sx)             # L̇ = (ΔX·ΔẊ)/L
    J .+= Diagonal(-2.0 .* N ./ L) * dL                 # length term
    J .+= Diagonal(p.Evec .* s ./ (L .^ 2)) * p.Sa      # area term (∂EA/∂x)
    return J
end

"""
    axial_stress_jacobian(t::DesignTangents, p::OptParams) -> Matrix

Exact Jacobian `∂σ/∂x` (n_elements × n_vars) of the per-element axial
stresses `σ = N/A`.
"""
function axial_stress_jacobian(t::DesignTangents, p::OptParams)
    res = t.res
    N = axial_force(res, p)
    J = axial_force_jacobian(t, p)
    J = Diagonal(1.0 ./ res.A) * J
    J .-= Diagonal(N ./ (res.A .^ 2)) * p.Sa
    return J
end
