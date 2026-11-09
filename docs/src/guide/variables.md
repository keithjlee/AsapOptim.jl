# Design variables

All variables share the pattern `Variable(target, start_value, lower_bound, upper_bound, ...)`. They have one of two semantics:

| Variable | Controls | Semantics |
|---|---|---|
| [`SpatialVariable`](@ref)`(node, value, lb, ub, axis)` | node position along `:X`/`:Y`/`:Z` | **additive**: the design entry is an *offset* from the modeled position (start at `0.0` to begin from the current geometry) |
| `SpatialVariable(node, vector, value, lb, ub)` | node position along an arbitrary directional "rail" | **additive**: `vector` is normalized, so the design entry is arc length along the rail |
| [`SectionVariable`](@ref)`(element, value, lb, ub, prop)` | a section property: `:A`, `:Ix`, `:Iy` or `:J` (flexural and torsional properties require a `FrameElement`) | **absolute**: the design entry *replaces* the property; the rest of the section is kept |
| [`AreaVariable`](@ref)`(element, value, lb, ub)` | cross-section area | **absolute**: the design entry *replaces* the area; `Ix`, `Iy`, `J` and the material are kept |
| [`JointVariable`](@ref)`(element, position, value, lb, ub)` | rotational end-spring stiffness `ky = kz` at `:start`, `:end` or `:both` of a `FrameElement` | **absolute** [force·length/rad]: semi-rigid connection design |
| [`QVariable`](@ref)`(element, value, lb, ub)` | force density of an `FDMelement` | **absolute**: for [force-density networks](networks.md) |

Notes:

- `AreaVariable` requires a geometric `Section`. A `RigiditySection` has no area to vary; parameterize its rigidities directly instead.
- Members may carry element loads under any variable: the fixed-end forces follow the design (see [Element loads](element-loads.md)).

## Coupling variables: symmetry and grouping

[`CoupledVariable`](@ref)`(target, parent, factor = 1.0)` makes `target` share `parent`'s design-vector entry, scaled by `factor`. This expresses equality constraints for free: the design space shrinks instead of gaining constraint rows.

```julia
# mirror-symmetric geometry: the right node's x-offset is the negative of the left's
x_left = SpatialVariable(left_node, 0.0, -1.0, 1.0, :X)
x_right = CoupledVariable(right_node, x_left, -1.0)

# member groups: all web elements share one area
a_web = AreaVariable(web_elements[1], 0.01, 1e-4, 0.05)
groups = [CoupledVariable(el, a_web) for el in web_elements[2:end]]

# joint pairs: couple another element's end to a JointVariable
# (the target may be an element, or an (element, position) tuple)
j = JointVariable(beam1, :start, 1e8, 1e5, 1e12)
mirror = CoupledVariable((beam2, :end), j)
```

Always couple to the *independent* parent; chains of couplings are rejected.

Variable types **mix freely in one problem**. Spatial, area, section-property and joint-stiffness variables, plus any couplings, can drive a single [`OptParams`](@ref) on models that mix frame and truss elements. Every combination is gradient-tested against finite differences in the test suite.
