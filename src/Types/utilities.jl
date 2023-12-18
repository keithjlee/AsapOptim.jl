"""
Set all values of a sparse matrix to *explicit* zeros. IE if S::SparseMatrixCSC has n non-zero values, explicit_zero(S) results in the same sparsity pattern with n 0s.
"""
explicit_zero(A::SparseMatrixCSC) = SparseMatrixCSC(size(A)..., A.colptr, A.rowval, zero(A.nzval))

function explicit_zero!(A::SparseMatrixCSC)
    A.nzval .= 0
end

"""
    get_inz(S::SparseMatrix, ids::Vector{Int64}, ndofpernode::Int64)
Get the indices in S.nzval, S.rowval with respect to a given set of indices in the global DOF order.
"""
function get_inz(S::SparseMatrixCSC{Float64, Int64}, ids::Vector{Int64}, ndofpernode::Int64)
    inzvals = Vector{Int64}()

    # the row/column index of the starting node DOF1 in S
    in1 = ids[1]
    # the row/column index of the ending node DOF2 in S
    in2 = ids[1+ndofpernode]

    offset = collect(-1:ndofpernode-2)

    for i in ids
        # the indices of S.nzval that correspond to non-zero values in column i
        rowrange = S.colptr[i]:S.colptr[i+1] - 1

        # the row indices of the non-zero values in S with respect to column i
        rvs = S.rowval[rowrange]

        # the index of node1 DOFS in S.nzval/S.rowval
        irv1 = findfirst(rvs .== in1) + rowrange[1] .+ offset

        # the index of node2 DOFS in S.nzval/S.rowval
        irv2 = findfirst(rvs .== in2) + rowrange[1] .+ offset

        inzvals = [inzvals; irv1; irv2]
    end

    inzvals
end

"""
    all_inz(model::Model)

For each element in the model, extract the location in the global stiffness matrix CSC structure corresponding to the column-wise values of the elemental stiffness matrix in GCS. IE given the result inz ∈ all_inz(model) that corresponds to element E:

```julia-repl   
julia> model.S.nzval[inz] = vec(E.K) + vec(E2.K) + ...
```

Where E2... are other elements that share nodes with E (may not exist).

This allows for a highly efficient method of regenerating the global stiffness matrix with new values while reusing the known sparsity pattern.
"""
all_inz(model::TrussModel) = [get_inz(model.S, id, 3) for id in getproperty.(model.elements, :globalID)]
all_inz(model::Model) = [get_inz(model.S, id, 6) for id in getproperty.(model.elements, :globalID)]

"""
get_local_global_DOF_activity(model::Asap.AbstractModel)
"""
function get_local_global_DOF_activity(model::Asap.AbstractModel)

    #indices of free/fixed DOFs
    i_free = model.freeDOFs
    i_fixed = model.fixedDOFs

    #get dictionary between freeDOFs and index in reduced K
    i_reduced = collect(1:length(i_free))
    init2reduced = Dict(i_free .=> i_reduced)

    #collectors
    local_ids = Vector{Vector{Int64}}()
    global_ids = Vector{Vector{Int64}}()

    for i in eachindex(model.elements)

        #current element
        element = model.elements[i]

        #element stiffness matrix in GCS
        k = element.K

        #associated global DOF indices
        dof_ids = element.globalID

        #which row/columns in k are active?
        i_active_local = [index for index in axes(k, 1) if !in(dof_ids[index], i_fixed)]

        #store
        push!(local_ids, i_active_local)

        #which DOFs in the global numbering system are active?
        i_active_dofs = dof_ids[i_active_local]

        #what indices of the reduced stiffness matrix are these associated with?
        i_active_dofs_reduced = [init2reduced[index] for index in i_active_dofs]

        #store
        push!(global_ids, i_active_dofs_reduced)

    end

    return local_ids, global_ids
end

"""
    get_inz_reduced(colptr::Vector{Int64}, rowval::Vector{Int64}, ids::Vector{Int64})

Get the indices in S.nzval, S.rowval with respect to a given set of active indices of an element.
"""
function get_inz_reduced(colptr::Vector{Int64}, rowval::Vector{Int64}, ids::Vector{Int64})
    inzvals = Vector{Int64}()

    for i in ids
        # the indices of K.nzval that correspond to non-zero values in column i
        rowrange = colptr[i]:colptr[i+1] - 1

        # the row indices of the non-zero values in S with respect to column i
        rvs = rowval[rowrange]

        # the index of the indices in ids in S.nzval/S.rowval
        idset = [findfirst(rvs .== val) for val in ids] .+ first(rowrange) .- 1

        inzvals = [inzvals; idset]
    end

    return inzvals
end
