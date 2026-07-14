#=
Gradient via ENZYME — LLVM-level AD, and (as of this writing) the FASTEST
backend for Asap: it beats even the legacy hand-differentiated
implementation on the publication benchmark problems.

Loading Enzyme alongside Asap activates `AsapEnzymeExt`, which imports the
two rules Enzyme cannot traverse (the sparse construction and the CHOLMOD
solve). Enzyme then needs TWO annotations, both shown below:

1. runtime activity — the objective captures the (constant) model/params,
   and Enzyme's static activity analysis cannot prove that's safe;
2. `function_annotation = Enzyme.Const` — tells Enzyme the objective
   closure's captured state is data, not something to differentiate.
   (Without it: `EnzymeMutabilityException`. Named top-level functions
   whose captured state is reached via globals also avoid the issue —
   `problem.jl` defines OBJ as a named function for exactly this reason.)

    ] add Enzyme DifferentiationInterface
=#

include("problem.jl")

using DifferentiationInterface
import Enzyme

backend = AutoEnzyme(;
    mode = Enzyme.set_runtime_activity(Enzyme.Reverse),
    function_annotation = Enzyme.Const,
)

prep = prepare_gradient(OBJ, backend, x0)

g = gradient(OBJ, prep, backend, x0)
println("∇C (Enzyme):   ", round.(g; sigdigits = 5))

t = time_gradient(() -> gradient(OBJ, prep, backend, x0))
println("median gradient time: ", round(t; digits = 3), " ms")
