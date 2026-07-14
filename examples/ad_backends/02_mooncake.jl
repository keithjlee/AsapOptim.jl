#=
Gradient via MOONCAKE — the actively developed successor to Zygote.

Loading Mooncake alongside Asap activates `AsapMooncakeExt`, which bridges
Asap's ChainRules rules into Mooncake's rule system (`@from_rrule`) and
teaches it that CHOLMOD factorizations are opaque solver state. After
that, nothing backend-specific is needed.

Mooncake compiles a rule per objective on first use (slower first call),
then typically runs with markedly LESS MEMORY than Zygote — and faster on
small problems.

    ] add Mooncake DifferentiationInterface
=#

include("problem.jl")

using DifferentiationInterface
import Mooncake

backend = AutoMooncake(; config = nothing)

# preparation does the rule compilation — do it ONCE per objective
prep = prepare_gradient(OBJ, backend, x0)

g = gradient(OBJ, prep, backend, x0)
println("∇C (Mooncake): ", round.(g; sigdigits = 5))

t = time_gradient(() -> gradient(OBJ, prep, backend, x0))
println("median gradient time: ", round(t; digits = 3), " ms")

#=
Gotcha worth knowing: if your objective closes over a model that was ALREADY
solved (its cache holds a CHOLMOD factorization), the extension's
`tangent_type` declaration for `CHOLMOD.Factor` is what keeps Mooncake from
trying to differentiate the factorization object itself.
=#
