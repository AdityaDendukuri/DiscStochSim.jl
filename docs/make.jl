using DiscStochSim
using Documenter

DocMeta.setdocmeta!(DiscStochSim, :DocTestSetup, :(using DiscStochSim); recursive=true)

makedocs(;
    modules=[DiscStochSim],
    authors="Aditya Dendukuri",
    sitename="DiscStochSim.jl",
    format=Documenter.HTML(;
        canonical="https://AdityaDendukuri.github.io/DiscStochSim.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AdityaDendukuri/DiscStochSim.jl",
    devbranch="main",
)
