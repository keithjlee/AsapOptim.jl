include("Utilities.jl")
export replace_values
export add_values

include("Geometry.jl")

include("Rtruss.jl")

include("Ktruss.jl")

include("K.jl")

include("Solve.jl")

include("Objective.jl")
export solve_truss
export compliance
export variation
export max_penalty
export min_penalty

export solve_network
export target

include("PostProcessing.jl")
export axial_force
export axial_stress