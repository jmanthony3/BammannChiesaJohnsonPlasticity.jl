using Documenter
using DocumenterCitations
using BammannChiesaJohnsonPlasticity
using Optimization, LossFunctions

DocMeta.setdocmeta!(BammannChiesaJohnsonPlasticity, :DocTestSetup, :(using BammannChiesaJohnsonPlasticity); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib"),
    style=:numeric
)

makedocs(;
    modules = [BammannChiesaJohnsonPlasticity, Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)],
    authors = "Joby M. Anthony III",
    repo    = "https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl/blob/{commit}{path}#{line}",
    sitename= "BammannChiesaJohnsonPlasticity.jl",
    doctest = false,
    format  = Documenter.HTML(;
        prettyurls  = get(ENV, "CI", "false") == "true",
        canonical   = "https://jmanthony3.github.io/BammannChiesaJohnsonPlasticity.jl",
        edit_link   = "main",
        assets      = String[],
    ),
    pages   = [
        "Home" => "index.md",
        "Bammann-Chiesa-Johnson Plasticity" => [
            "Base Package" => "BammannChiesaJohnsonPlasticity.md",
            "Metals" => [
                "Metal Specific Types" => "Metals.md",
                "Bammann1990Modeling" => "Bammann1990Modeling.md",
            ]
        ],
        "Extensions" => [
            "Extending Functionality" => "extensions.md",
            "Optimization.jl" => "ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt.md"
        ]
    ],
    plugins = [bib],
)

deploydocs(;
    repo="github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl",
    devbranch="main",
)
