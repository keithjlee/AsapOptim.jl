# full stiffness matrix
K = copy(prob.K)

# vector of sorted active DOF indices
inds = copy(params.freeids)

inds_fixed = model.fixedDOFs

#reduced indices
ired = collect(1:length(inds))

#ind comparison
icomp = [inds ired]

# dictionary
init2reduced = Dict(inds .=> ired)

# target reduced stiffness matrix
Kred = copy(K[inds, inds])

#=
 indices of points in inds where inds[i+1] - inds[i] > 1.

 IE the DOFs between 
=#
i_startof_skip = findall(diff(inds) .!= 1)

#=

INEFFICIENT TEST 1: POPULATE A NEW REDUCED STIFFNESS MATRIX USING THE INDEX DICTIONARY

=#

# approach 1: append to initialized sparse matrix K
K1 = spzeros(Float64, length(ired), length(ired))

local_ids = Vector{Vector{Int64}}()
global_ids = Vector{Vector{Int64}}()

for i = 1:model.nElements

    # stiffness matrix in GCS for element i
    k = model.elements[i].K

    # DOF of element ends [i_start_x, i_start_y, i_start_z, i_end_x, i_end_y, i_end_z]
    dofids = model.elements[i].globalID

    # startindex = dofids[1]
    # endindex = dofids[4]

    # active_local = [
    #     in(startindex, inds_fixed) ? [] : 1:3;
    #     in(endindex, inds_fixed) ? [] : 4:6
    # ]

    active_local = [index for index in 1:6 if !in(dofids[index], inds_fixed)]
    
    push!(local_ids, active_local)

    active_dofs = dofids[active_local]
    active_dofs_reduced = [init2reduced[i] for i in active_dofs]

    push!(global_ids, active_dofs_reduced)
    #directly append to Knew ... removes 0s by default
    @views K1[active_dofs_reduced, active_dofs_reduced] .+= k[active_local, active_local]
end




iel = rand(1:model.nElements)

gid = global_ids[iel]
lid = local_ids[iel]

i = 2
gid_single = gid[i]

rowrange = Kred.colptr[gid_single]:Kred.colptr[gid_single+1] - 1
rowvals = Kred.rowval[rowrange]

idset = [findfirst(rowvals .== val) for val in gid]

@time inz_ids = [get_inz_reduced(Kred, id) for id in global_ids]

K2 = copy(Kred)
fill!(nonzeros(K2), 0.)

for i in eachindex(model.elements)

    iloc = local_ids[i]

    K2.nzval[inz_ids[i]] .+= model.elements[i].K[iloc, iloc][:]

end