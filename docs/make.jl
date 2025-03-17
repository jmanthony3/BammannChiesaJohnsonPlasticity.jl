using Pkg; Pkg.precompile()
using Documenter
using DocumenterCitations
using BammannChiesaJohnsonPlasticity

DocMeta.setdocmeta!(BammannChiesaJohnsonPlasticity, :DocTestSetup, :(using BammannChiesaJohnsonPlasticity); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib"),
    style=:numeric
)

makedocs(;
    modules = [BammannChiesaJohnsonPlasticity],
    authors = "Joby M. Anthony III, Julian Tse Lop Kun",
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
        "Bammann-Chiesa-Johnson Plasticity" => "bcjplasticity.md",
        "Metals" => "metals.md",
        "Extensions" => "extensions.md"
    ],
    plugins = [bib],
)

deploydocs(;
    repo="github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl",
    devbranch="main",
)
