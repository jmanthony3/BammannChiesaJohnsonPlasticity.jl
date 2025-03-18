# Base

```@meta
CurrentModule = BammannChiesaJohnsonPlasticity
```

This package adds overloads and sub-types from the [`ContinuumMechanicsBase.jl`](https://github.com/TRACER-LULab/ContinuumMechanicsBase.jl.git) package for the BCJ-plasticity model.
Most importantly, the `AbstractBCJModel` and `AbstractBCJTest` sub-types are defined first such that subsequent types and functions might dispatch under these two types.
Some helper functions are also defined for calculating the mathematical norm from a vector representing the upper triangle of a symmetric tensor and the von Mises equivalent of a symmetric tensor.
The elements of this vector for the upper triangle of a symmetric vector are in row-major order.
That is, the upper triangle of a symmetric, second-rank tensor may be represented as a vector: e. g. $\underbar{\underbar{A}} \vec{\equiv} [A_{11}, A_{12}, A_{13}, A_{22}, A_{23}, A_{33}]$.
Currently, the package includes support for the following categories of materials:

```@contents
Modules = [BammannChiesaJohnsonPlasticity]
Pages   = [
    "Metals.md"]
Depth   = 1
```

## Types
```@autodocs
Modules = [BammannChiesaJohnsonPlasticity]
Order   = [:type]
Pages   = ["BammannChiesaJohnsonPlasticity.jl"]
```

## Functions
```@autodocs
Modules = [BammannChiesaJohnsonPlasticity]
Order   = [:function]
Pages   = ["BammannChiesaJohnsonPlasticity.jl"]
```

## Index
```@index
Modules = [BammannChiesaJohnsonPlasticity]
Order   = [:type, :function]
Pages   = ["BammannChiesaJohnsonPlasticity.md"]
```
