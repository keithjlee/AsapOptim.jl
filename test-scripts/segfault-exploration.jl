# activate environment and load dependencies
using AsapOptim
import AsapOptim: Asap

#=
CURRENT MWE FOR SEGFAULT

Repeatedly execute everything below starting at `n = ...` in the REPL (i.e. press shift + enter).
If you reach the end of the script, try changing the value of n to some even integer value (e.g. n = 26, 40, 30, 18, ...) and rerun the script.
Eventually, a fatal crash (segfault) will occur. This may happen on the first time you run this script, or the fifth, etc.

Notes:
- Generally, this occurs sooner with larger values of n.
- Occurs in all version of Julia 1.9 and higher (this environment was made in 1.10.3). DOES NOT OCCUR in the latest version of Julia 1.8 (1.8.5).
- Occurs consistently on a M1 Pro Macbook Pro and a Windows 11 Intel PC.
- This does NOT occur if you enclose everything in a for loop and run for a bunch of iterations
=#


# select some even integer value for n (number of nodes in the square grid)
n = 40

# [Asap] make a cross section 
section = Asap.Section(
    20e-3, # A [m²]
    2e8, # E [kN/m²]
    8e7, # G [kN/m²]
    7.95e-4, # Ix [m⁴]
    9.2e-5, # Iy [m⁴]
    3.11e-6 # J [m⁴]
)

# hyper parameters for grid
begin
    Lx = 25.
    Ly = 15.
end


# generate a square grid structural model
begin

    # applied load vector
    load = [0., 0., -20] #kN

    # nodal positions
    x_positions = range(0, Lx, n)
    y_positions = range(0, Ly, n)

    # distance between nodes
    dx = Lx / (n-1)
    dy = Ly / (n-1)

    # collectors
    xyz = Vector{Vector{Float64}}()
    Xmatrix = zeros(Float64, n, n)
    Ymatrix = zeros(Float64, n, n)
    igrid = zeros(Int64, n, n)

    index = 1
    for iy = 1:n
        for ix = 1:n

            x = x_positions[iy]
            y = y_positions[ix]

            igrid[ix, iy] = index
            index += 1

            push!(xyz, [x, y, 0.])
            Xmatrix[ix, iy] = x
            Ymatrix[ix, iy] = y

        end
    end

    # extract corner node indices to fix
    support_indices = [igrid[1, 1], igrid[n, 1], igrid[1, n], igrid[n, n]]

    # [Asap] make nodes
    nodes = [Asap.Node(pos, :free, :free) for pos in xyz]

    # [Asap] make support nodes
    for node in nodes[support_indices]
        Asap.fixnode!(node, :pinned)
        node.id = :support
    end

    # [Asap] element collector
    elements = Vector{Asap.Element}()

    #horizontal elements
    for i = 1:n
        for j = 1:n-1
            index = [igrid[i,j], igrid[i,j+1]]
            push!(elements, Asap.Element(nodes[index]..., section))
        end
    end

    #vertical elements
    for j = 1:n
        for i = 1:n-1
            index = [igrid[i,j], igrid[i+1,j]]
            push!(elements, Asap.Element(nodes[index]..., section)) 
        end
    end

    # [Asap] loads
    loads = [Asap.NodeForce(node, load) for node in nodes[:free]]

    # [Asap] assemble model and solve
    model = Asap.Model(nodes, elements, loads)
    Asap.solve!(model)
end;

#=
design variable indices

This does not make much sense in this MWE, but is useful when using CoupledVariables to enforce symmetry of nodal positions
=#
begin
    @assert n % 2 == 0

    igrid = reshape(1:n^2, n, n)

    imid = Int(n / 2)

    iparent = igrid[2:imid, 2:imid]

    ichild1 = reverse(igrid[2:imid, imid+1:end-1], dims = 2)
    factors1 = [-1., 1.]

    ichild2 = reverse(igrid[imid+1:end-1, 2:imid], dims = 1)
    factors2 = [1., -1.]

    ichild3 = reverse(igrid[imid+1:end-1, imid+1:end-1])
    factors3 = [-1., -1.]
end

# [AsapOptim] make design variables
begin
    vars = Vector{FrameVariable}()

    fac = .9
    x = dx * fac / 2
    y = dy * fac / 2
    z = 1.5


    for i in eachindex(iparent)

        i0 = iparent[i]
        i1 = ichild1[i]
        i2 = ichild2[i]
        i3 = ichild3[i]

        push!(vars, SpatialVariable(i0, 0., -x, x, :X))
        push!(vars, SpatialVariable(i1, 0., -x, x, :X))
        push!(vars, SpatialVariable(i2, 0., -x, x, :X))
        push!(vars, SpatialVariable(i3, 0., -x, x, :X))

        push!(vars, SpatialVariable(i0, 0., -y, y, :Y))
        push!(vars, SpatialVariable(i1, 0., -y, y, :Y))
        push!(vars, SpatialVariable(i2, 0., -y, y, :Y))
        push!(vars, SpatialVariable(i3, 0., -y, y, :Y))
    end
end

# [AsapOptim] make optimization parameters
FrameOptParams(model, vars)

#=
If you've come this far without a segfault, congrats! Trying changing the value of `n` above or changing some hyperparameters (Lx, x, y, ....) and rerunning.
=#