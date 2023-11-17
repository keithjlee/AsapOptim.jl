## On general displacement adjoint
u = K⁻¹P

if obj = f(u), then the gradient of obj with respect to an independent variable x is achieved through the chain rule:

dObj/dx = ... ⋅ du/dK ⋅ df/du = ... ⋅ du/dK ⋅ ū

For this rule, we are concerned with finding du/dK, or the [ndof × ndof] matrix of sensitivites that we can propagate backwards to the final objective.

Given df/du = [ndof × 1] = ū is the gradient of the objective function with respect to displacements u, the sensitivity is:

du/dK = - uᵀ ⊗ K⁻¹
df/dK = du/dK ū = - (uᵀ ⊗ K⁻¹)ū

Can be rearranged such that:
ΔK = K⁻¹ū

df/dK = -uᵀ ⊗ ΔK

Which is an [ndof × ndof] matrix where:

Columnᵢ = uᵢ .* ΔK