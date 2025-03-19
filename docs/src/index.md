# BammannChiesaJohnsonPlasticity
Documentation for [BammannChiesaJohnsonPlasticity](https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl).

## Purpose
The `BammannChiesaJohnsonPlasticity.jl` package is a Julian implementation of Bammann-Chiesa-Johnson (BCJ) plasticity from select publications.
**The goal of this package is to provide a repository and interface for the BCJ-plasticity model and its equations to implement and interact with those equations from publications, as well as, give opportunity for developing or improving the model.**
The BCJ model has been applied to many materials over the years; however, only those models from a select few papers for modeling the inelastic deformation of _metals_ are currently defined in this package.
It is a hope of the authors that published models expanding, refining, or applying the BCJ model onto a wider selection of materials might be included in this package along with the model constants employed in those publications for the appropriate chemical system.
Having such a repository for versions of the BCJ-plasticity model, variations of coupling BCJ-plasticity with other ISV models (damage, recrystallization, etcetera), and the relevant constants maintains transparency and reproducibility of publications.
Many publications discuss the physics and mechanics of the BCJ-plasticity model, but a small primer is included in the next section for those unfamiliar with the model.
At this time, all other sections after the primer will only include comments on the behavior of the package and not the limitations of the model as this is discussed in the literature.

## Methodology
This `BCJPlasticity.jl` package employs and overloads the [`ContinuumMechanicsBase.jl`](https://github.com/TRACER-LULab/ContinuumMechanicsBase.jl.git) package for the types and functions common to continuum mechanics modeling.
Currently, the package acts as a point simulator to update the material states and assumes a Poisson's ratio of 0.5.
Although the package can be used out-of-the-box to plot the specified BCJ model with the appropriate loading conditions, material properties, and model constants, the package can be extended to calibrate the constants against experimental data via the [`Optimization.jl`](https://github.com/SciML/Optimization.jl.git) package.
For an example of interacting with and calibrating BCJ model constants, the reader is referred to the notebooks folder for a set of [`Pluto.jl`](https://github.com/fonsp/Pluto.jl.git) notebooks as examples.
Furthermore, in that folder contains a notebook for interacting with the [Johnson and Cook (1983)](@cite johnsonCONSTITUTIVEMODELDATA1983a).
This `JohnsonCook.jl` notebook `include`s a neighboring script (`JohnsonCook-functions.jl`) which defines all the necessary types and functions if one wants to create or modify a viscoplasticity model for interaction or calibration as in a Pluto notebook.
The empirical Johnson-Cook model is included as an example given that its widespread recognition might minimize the learning curve for encoding a material model into a Julian implementation appropriate for calibration with `Optimization.jl` or interaction in a Pluto notebook.

## Nomenclature
As much as possible, equations listed in the documentation for any version of the BCJ-plasticity model will attempt to match the form of the equations from the cited publication.
For consistency, some cosmetics are applied in the documentation to maintain clarity as to the rank of tensors.
For English letters, first and second rank tensors will be denoted as lower and uppercase, bold-face symbols, respectively: e. g. ``\mathbf{x}`` is a first rank tensor while ``mathbf{A}`` is second rank.
For Greek letters, Einsteinian notation, instead of bold type facing, is used to connote tensors: e. g. ``\underset{\sim}{\sigma}`` could be the second rank tensor for stress.
Many functions used throughout this package, as well as its overloads for `ContinuumMechanicsBase.jl`, assign the prescribed BCJ-plasticity model to the variable `ψ` and the applied boundary conditions assigned to `Ω`.

## References
```@bibliography
```
