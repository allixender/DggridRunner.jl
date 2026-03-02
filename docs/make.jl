using Documenter
using Literate
using DGGRIDRunner
import DGGRIDRunner.DGGRIDParams
import DGGRIDRunner.AuthalicConversion

# ---------------------------------------------------------------------------
# Generate example pages from Literate.jl source files
# ---------------------------------------------------------------------------

examples_dir = joinpath(@__DIR__, "..", "examples")
examples_out = joinpath(@__DIR__, "src", "examples")
mkpath(examples_out)

for example in [
    "whole_earth_igeo7.jl",
    "coarse_cells_workflow.jl",
    "clip_cells_workflow.jl",
]
    Literate.markdown(
        joinpath(examples_dir, example), examples_out;
        execute   = false,
        codefence = "```julia" => "```",
    )
end

# Copy repo-root assets referenced from docs/src/index.md
let src = joinpath(@__DIR__, "..", "day-04-hexa.png")
    isfile(src) && cp(src, joinpath(@__DIR__, "src", "day-04-hexa.png"); force = true)
end

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

makedocs(
    sitename = "DGGRIDRunner.jl",
    authors  = "Alexander Kmoch and contributors",
    modules  = [DGGRIDRunner, DGGRIDParams, AuthalicConversion],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home"          => "index.md",
        "Examples" => [
            "examples/whole_earth_igeo7.md",
            "examples/coarse_cells_workflow.md",
            "examples/clip_cells_workflow.md",
        ],
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/allixender/DggridRunner.jl.git",
    devbranch = "main",
    push_preview = false,
)
