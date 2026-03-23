using Documenter, BroadbandBPM

makedocs(
    modules = [BroadbandBPM],
    pages = [
        "Home" => "index.md",
    ],
    sitename="BroadbandBPM.jl",
    authors="Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/BroadbandBPM.jl.git"
)