#=
Gradient via ZYGOTE — the incumbent.

Zygote needs nothing special: loading it activates Asap's ChainRules
extension, and every rule the pipeline needs is already registered.

    ] add Zygote DifferentiationInterface
=#

include("problem.jl")

using DifferentiationInterface
import Zygote

backend = AutoZygote()

# one-off gradient
g = gradient(OBJ, backend, x0)
println("∇C (Zygote):   ", round.(g; sigdigits = 5))

# for repeated evaluation (optimization loops), prepare once:
prep = prepare_gradient(OBJ, backend, x0)
t = time_gradient(() -> gradient(OBJ, prep, backend, x0))
println("median gradient time: ", round(t; digits = 3), " ms")

# Zygote also works directly, without DifferentiationInterface:
g_direct = Zygote.gradient(OBJ, x0)[1]
@assert g ≈ g_direct
