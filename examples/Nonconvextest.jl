using Asap, AsapToolkit, AsapOptim
using kjlMakie; set_theme!(kjl_dark)
import Nonconvex
Nonconvex.@load NLopt


### Create a spaceframe
#meta parameters
begin
    nx = 21
    dx = 1200.
    ny = 15
    dy = 1000.
    dz = 2000.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

# generate and extract model
sf = generatespaceframe(nx, 
    dx, 
    ny, 
    dy, 
    dz,
    tube,
    ; 
    load = [0., 0., -30e3], 
    support = :y);

model = sf.truss;

for node in model.nodes[sf.iY1]
    fixnode!(node, :zfixed)
end
updateDOF!(model); solve!(model)

begin
    dfac = Observable(1.)
    p0 = @lift(Point3.(getproperty.(model.nodes, :position) .+ $dfac * getproperty.(model.nodes, :displacement)))
    e0 = @lift(vcat([$p0[id] for id in getproperty.(model.elements, :nodeIDs)]...))
    f0 = getindex.(getproperty.(model.elements, :forces), 2)
    c0 = maximum(abs.(f0)) .* (-1, 1) .* .2
    s0 = @lift($p0[findall(model.nodes, :support)])

    fig = Figure(
        # backgroundcolor = :blue
        )
    ax = Axis3(fig[1,1],
        protrusions = 75,
        aspect = :data)

    gridtoggle!(ax); simplifyspines!(ax)

    linesegments!(e0,
        colormap = pink2blue,
        color = f0,
        colorrange = c0,
        linewidth = 4)

    scatter!(s0,
        color = :yellow,
        strokecolor = :black,
        strokewidth = 4,
        size = 100)
    

    fig
end

#make variables
begin
    vars = Vector{TrussVariable}()

    # nodal Z variables
    lb = -500.
    ub = 3500.

    lbb = -2000.
    ubb = 750.

    #nodal xy variables for bottom nodes
    lbxy = -750.
    ubxy = 750.

    #spatial variables
    for node in model.nodes

        #top nodes can move in X, Y, Z
        if node.id == :top
            push!(vars, SpatialVariable(node, 0., lb, ub, :Z))
            # push!(vars, SpatialVariable(node, 0., lbxy, ubxy, :X))
            # push!(vars, SpatialVariable(node, 0.,  lbxy, ubxy, :Y))
        end

        #bottom nodes can move in Z
        if node.id == :bottom
            push!(vars, SpatialVariable(node, 0., lbb, ubb, :Z))
        end
    end

    for el in model.elements
        push!(vars, AreaVariable(el, 15_000., 1_000., 30000.))
    end

end

# generate a problem
params = TrussOptParams(model, vars);

# define objective function
function obj(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = addvalues(p.X, p.indexer.iX, values[p.indexer.iXg])
    Y = addvalues(p.Y, p.indexer.iY, values[p.indexer.iYg])
    Z = addvalues(p.Z, p.indexer.iZ, values[p.indexer.iZg])
    A = replacevalues(p.A, p.indexer.iA, values[p.indexer.iAg])

    # vₑ
    v = AsapOptim.getevecs(X, Y, Z, p)

    # Lₑ
    l = AsapOptim.getlengths(v)

    dot(A, l)
end

function cstr(values::Vector{Float64}, p::TrussOptParams)

    res = solvetruss(values, p)

    stress = maximum(axialstress(res, p)) - 350
    disp = - nx * dx / 360 - minimum(res.U[3:3:end])

    (stress, disp)
    # stress
end

#test
@time obj(params.values, params)

@time cstr(params.values, params)

F = Nonconvex.TraceFunction(x -> obj(x, params))
optmodel = Nonconvex.Model(F)

#add variables and constraints
Nonconvex.addvar!(optmodel, params.lb, params.ub, init = copy(params.values))
Nonconvex.add_ineq_constraint!(optmodel, x -> cstr(x, params))

#define algorithm
@time begin
    alg = NLoptAlg(:LD_MMA)
    opts = NLoptOptions(ftol_rel = 1e-6, maxeval = 500)
    res = Nonconvex.optimize(optmodel, alg, params.values, options = opts)
end

@show res.minimum
@show cstr(res.minimizer, params)

m2 = updatemodel(params, res.minimizer);

#solution
begin
    sol = res.minimizer

    a = AsapOptim.replacevalues(params.A, params.indexer.iA, sol[params.indexer.iAg])

    lfac = Observable(5.)
    lw = @lift(a ./ maximum(a) .* $lfac)
end

ms = 20
nb = 25
begin
    p1 = @lift(Point3.(getproperty.(m2.nodes, :position) .+ $dfac * getproperty.(m2.nodes, :displacement)))
    e1 = @lift(vcat([$p1[id] for id in getproperty.(m2.elements, :nodeIDs)]...))
    f1 = getindex.(getproperty.(m2.elements, :forces), 2)
    c1 = maximum(abs.(f1)) .* (-1, 1) .* .2
    s1 = @lift($p1[findall(m2.nodes, :support)])

    fig = Figure()

    ### geometry plots
    ax = Axis3(fig[1,1],
        aspect = :data)

    gridtoggle!(ax); simplifyspines!(ax)

    oldpos = linesegments!(e0,
        colormap = pink2blue,
        color = f0,
        colorrange = c0,
        linewidth = 4)

    scatter!(s0,
        color = :yellow,
        strokecolor = :black,
        strokewidth = 4,
        size = 100)


    ax2 = Axis3(fig[1,2],
        aspect = :data)

    hidedecorations!(ax2); hidespines!(ax2)

    newpos = linesegments!(e1,
        color = f1,
        colormap = pink2blue,
        colorrange = c0,
        linewidth = lw)
    
    linkaxes!(ax, ax2)

    ### AL scatters
    ax_al = Axis(fig[2,1],
        aspect = 4,
        xlabel = "A [mm²]",
        ylabel = "L [mm]")

    scatter!(getproperty.(getproperty.(model.elements, :section), :A),
        getproperty.(model.elements, :length),
        color = sign.(f0),
        colormap = pink2blue,
        markersize = ms)

    ax2_al = Axis(fig[2,2],
        aspect = 4,
        xlabel = "A [mm²]",
        ylabel = "L [mm]")

    scatter!(a,
        getproperty.(m2.elements, :length),
        color = sign.(f1),
        colormap = pink2blue,
        markersize = ms)

    linkyaxes!(ax_al, ax2_al)

    ### stress histograms
    ax_stress = Axis(fig[3,1],
        aspect = 3,
        xlabel = "σ [MPa]",
        )

    hideydecorations!(ax_stress)

    hist!(f0 ./ sec.A,
        color = :white)

    ax2_stress = Axis(fig[3,2],
        aspect = 3,
        xlabel = "σ [MPa]",
        )

    hideydecorations!(ax2_stress)

    hist!(f1 ./ a,
        color = :white)

    on(dfac) do _
        reset_limits!(ax)
        reset_limits!(ax2)
    end

    fig
end
