# Extensions
The functionality of the base package can be extended upon loading the loading other packages.
Currently, the following extensions can be loaded for additional interaction with the BCJ-plasticity model, and the reader is directed to those pages for the extension of interest.
Future work may include an extension for the [`Ferrite.jl`](https://github.com/Ferrite-FEM/Ferrite.jl.git) package to perform finite element analyses in pure Julia with the BCJ-plasticity model.

```@contents
Modules = [
    Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)]
Pages   = [
    "ext/OptimizationBCJPlasticityExt.md"]
Depth   = 1
```
