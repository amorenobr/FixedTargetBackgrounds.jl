using FixedTargetBackgrounds
using Documenter

DocMeta.setdocmeta!(FixedTargetBackgrounds, :DocTestSetup, :(using FixedTargetBackgrounds); recursive=true)

makedocs(;
    modules=[FixedTargetBackgrounds],
    authors="Alexander Moreno Briceño <alexander.moreno@uan.edu.co> and contributors",
    sitename="FixedTargetBackgrounds.jl",
    format=Documenter.HTML(;
        canonical="https://amorenobr.github.io/FixedTargetBackgrounds.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/amorenobr/FixedTargetBackgrounds.jl",
    devbranch="main",
)
