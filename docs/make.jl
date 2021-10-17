push!(LOAD_PATH, joinpath("..", "src"))

using Documenter, FSPS

makedocs(
    sitename = "FSPS.jl",
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/elijahmathews/FSPS.jl.git",
    devbranch = "primary",
)
