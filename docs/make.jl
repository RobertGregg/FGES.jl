using FGES
using Documenter

DocMeta.setdocmeta!(FGES, :DocTestSetup, :(using FGES); recursive=true)

makedocs(;
    modules=[FGES],
    authors="RobertGregg <robwgregg@gmail.com> and contributors",
    repo="https://github.com/RobertGregg/FGES.jl/blob/{commit}{path}#{line}",
    sitename="FGES.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
