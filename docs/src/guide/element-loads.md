# Element loads

Members may carry element loads: `LineLoad`/`DistributedLoad`, `PointLoad`, `PointMoment` and `SelfWeight`. Asap lowers them to fixed-end forces, and **those forces depend on the design**:

| Fixed-end forces depend on | Through | Variables that reach them |
|---|---|---|
| member length and orientation | integration along the member; local–global rotation | [`SpatialVariable`](@ref) on the member's nodes |
| ρA | `SelfWeight` | [`AreaVariable`](@ref), `SectionVariable(:A)` |
| end springs and EI | condensation of released or semi-rigid ends | [`JointVariable`](@ref), `SectionVariable(:Ix/:Iy)` |

Since AsapOptim 1.2 (with Asap 1.2), the load vector is re-evaluated at every design, `F(x) = P − Pf(x)`. Gradients therefore include the load path, and they match finite differences of rebuilt, re-solved models. Earlier versions froze the load vector at the base design, which gave silently wrong gradients for all three rows above.

## Semantics under design change

- `LineLoad` / `DistributedLoad` are **per-unit-length intensities** over fractions of the member, so the total load follows the member's current length. This suits self-weight-like and tributary loads.
- `PointLoad` / `PointMoment` keep their **fractional position**, so they ride with the member.
- `SelfWeight` follows the current ρA and length.
- `:local` load directions rotate with the member; `:global` ones do not.
- Nodal loads (`NodeForce`, `NodeMoment`) are constant.

The condensed end moments of a uniformly loaded beam with symmetric end springs κ are `(wL²/12) · κL / (κL + 2EI)`. They run from 0 (hinged) to `wL²/12` (clamped) across the range a joint optimizer explores, so joint-stiffness design on gravity-loaded frames needs this load path.

## Member forces

`axial_force(res, params)` adds the element loads' axial fixed-end force to the elastic term. This matches Asap's stateful `axial_force(res, el)`: the axial action at the end node, tension-positive. It matters for inclined loaded members and for columns under self-weight. `axial_stress` and the implicit Jacobians ([`axial_force_jacobian`](@ref), [`axial_stress_jacobian`](@ref)) include the same term.

## Keeping the loads fixed

For studies where the load is fixed in space, compile the parameters with `loads = :frozen`:

```julia
params = OptParams(model, vars; loads = :frozen)
```

The base-design load vector `params.F` is then used for every design. Problems without element loads behave identically either way and pay nothing extra.
