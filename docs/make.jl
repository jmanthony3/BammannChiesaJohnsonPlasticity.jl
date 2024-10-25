using Pkg; Pkg.precompile()
using Documenter
using BammannChiesaJohnsonPlasticity

DocMeta.setdocmeta!(BammannChiesaJohnsonPlasticity, :DocTestSetup, :(using BammannChiesaJohnsons); recursive=true)

makedocs(;
    modules=[BammannChiesaJohnsonPlasticity],
    authors="Joby M. Anthony III, Julian Tse Lop Kun",
    repo="https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl/blob/{commit}{path}#{line}",
    sitename="BammannChiesaJohnsonPlasticity.jl",
    doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmanthony3.github.io/BammannChiesaJohnsonPlasticity.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl",
    devbranch="main",
)
