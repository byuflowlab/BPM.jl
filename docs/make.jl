using Documenter, BPM

makedocs(
    modules = [BPM],
    pages = [
        "Home" => "index.md",
    ],
    sitename="BPM.jl",
    authors="Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/BPM.jl.git"
)