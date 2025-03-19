# OptimizationBCJPlasticityExt

```@meta
CurrentModule = Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)
```

This extension is meant to provide additional functionality to the base package with the intent of performing optimization studies with the [`Optimization.jl`](https://github.com/SciML/Optimization.jl.git) and [`LossFunctions.jl`](https://github.com/JuliaML/LossFunctions.jl.git) packages to calibrate material model constants against experimental data.
This interface with `Optimization.jl` is largely accomplished by overloading three functions from `ContinuumMechanicsBase.jl` with the newly defined types in the `BCJPlasticity.jl` package:

- `parameters`,
- `parameter_bounds`, and
- `MaterialOptimizationProblem`

Most significant of these overloads is that for `parameter_bounds` wherein the default value for `lb` of the returned tuple, which represents the lower bounds of possible values during optimization, is constrained to be ``[0, \infty)`` to maintain physical admissibility and self-consistency in the constitutive and ISV rate equations.