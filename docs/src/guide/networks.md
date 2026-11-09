# Force-density networks

The same variables → parameters → solve pattern works on `Asap.Network` (force-density form-finding). [`solve_network`](@ref) is differentiable in the force densities end to end:

```julia
qv = QVariable(network.elements[1], 1.0, 0.1, 10.0)
vars = AbstractVariable[qv; [CoupledVariable(el, qv) for el in network.elements[2:5]]]

nparams = NetworkOptParams(network, vars)

function obj(x)
    res = solve_network(x, nparams)          # NetworkOptResults: X, Y, Z, Q, L
    sum(abs2, member_forces(res))            # forces are Q .* L, tension-positive
end

Zygote.gradient(obj, copy(nparams.values))
```

[`NetworkOptParams`](@ref) rejects mixed (per-axis) node fixity: the differentiable force-density path solves all three coordinates against one partition.
