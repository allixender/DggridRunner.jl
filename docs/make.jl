using Documenter
using DGGRIDRunner
import DGGRIDRunner.DGGRIDParams
import DGGRIDRunner.AuthalicConversion

makedocs(
    sitename = "DGGRIDRunner.jl",
    authors = "Alexander Kmoch and contributors",
    modules = [DGGRIDRunner, DGGRIDParams, AuthalicConversion],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home"          => "index.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/allixender/DggridRunner.jl.git",
    devbranch = "main",
    push_preview = false,
)
