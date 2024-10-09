using DUNEatLapp
using Documenter

DocMeta.setdocmeta!(DUNEatLapp, :DocTestSetup, :(using DUNEatLapp); recursive=true)

makedocs(;
    modules=[DUNEatLapp],
    authors="Mael <mael.martin@lapp.in2p3.fr>",
    sitename="DUNEatLapp.jl",
    format=Documenter.HTML(;
        canonical="https://MaelMartin17.github.io/DUNEatLapp.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MaelMartin17/DUNEatLapp.jl",
    devbranch="main",
)
