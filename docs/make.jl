# Build the AsapOptim.jl documentation site.
#
# Local build (from the repository root):
#   julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
#   julia --project=docs docs/make.jl
# then open docs/build/index.html (prettyurls are off outside CI, so direct
# file browsing works). In the AsapEcosystem workspace, also
# `Pkg.develop(PackageSpec(path="../Asap"))` to build against the sibling Asap.

using Documenter, AsapOptim, Asap
using Zygote   # activates Asap's ChainRules extension for the guide's @example gradients

makedocs(
    sitename = "AsapOptim.jl",
    authors = "Keith J. Lee",
    modules = [AsapOptim],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://keithjlee.com/AsapOptim.jl",
        edit_link = "main",
        size_threshold = 400 * 2^10,
        size_threshold_warn = 250 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "Quick start" => "guide/quickstart.md",
            "Design variables" => "guide/variables.md",
            "Evaluating designs" => "guide/evaluating.md",
            "Element loads" => "guide/element-loads.md",
            "Constrained optimization and AD" => "guide/constrained.md",
            "Force-density networks" => "guide/networks.md",
            "Examples" => "guide/examples.md",
        ],
        "Migrating from v0.1" => "migration.md",
        "API reference" => [
            "Variables" => "api/variables.md",
            "Parameters" => "api/parameters.md",
            "Results and evaluation" => "api/results.md",
            "Implicit differentiation" => "api/implicit.md",
            "Force-density networks" => "api/networks.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/keithjlee/AsapOptim.jl",
    devbranch = "main",
    push_preview = false,
)
