# BammannChiesaJohnsonPlasticity
Documentation for [BammannChiesaJohnsonPlasticity](https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl).

The `BammannChiesaJohnsonPlasticity.jl` package is a Julian implementation of Bammann-Chiesa-Johnson (BCJ) plasticity to model the inelastic deformation of solid materials at a given temperature and strain rate with Internal State Variables (ISVs).
The BCJ-plasticity model satisfies the _Clausius-Duhem inequality_ by placing thermodynamic constraints on ISVs and their rate equations ([Coleman and Gurtin (1967)]([colemanThermodynamicsInternalState1967](@cite))) which are motivated by equations for dislocation mechanics.
The model assumes that for any continuum mechanics configuration each material point at each moment in time satisfies some equilibrium criterion for a triplet point composed of the deformation gradient, temperature, and current set of ISVs.
The ISVs in the current configuration intrinsically represent the evolution of the ISVs because of the material history afforded by the rate equations.
These ISVs and their rate equations are used to increment the current configuration from the previous configuration by a trial guess that the material deforms elastically and correcting this trial guess for any plastic deformation accrued in the increment: this is the scheme of the _Radial-Return Method_ from [Krieg and Krieg (1977)]([kriegAccuraciesNumericalSolution1977](@cite)).

**The goal of this package is to provide an interface for the BCJ-plasticity model and its equations to implement and interact with those equations from publications, as well as, give opportunity for developing or improving the model.**
The BCJ model has been applied to many materials over the years; however, only those models from a select few papers for modeling the inelastic deformation of _metals_ are currently defined in this package.
It is a hope of the authors that published models expanding, refining, or applying the BCJ model onto a wider selection of materials might be included in this package along with the model constants employed in those publications for the appropriate chemical system.

This package employs the [`ContinuumMechanicsBase.jl`](https://github.com/TRACER-LULab/ContinuumMechanicsBase.jl.git) package for the types and functions common to continuum mechanics modeling.
Currently, the package acts as a point simulator to update the material states and assumes a Poisson's ratio of 0.5.
Although the package can be used out-of-the-box to plot the specified BCJ model with the appropriate loading conditions, material properties, and model constants, the package can be extended to calibrate the constants against experimental data via the [`Optimization.jl`](https://github.com/SciML/Optimization.jl.git) package.
For an example of interacting with and calibrating BCJ model constants, the reader is referred to the [notebooks](../../notebooks/) folder for a set of [`Pluto.jl`](https://github.com/fonsp/Pluto.jl.git) notebooks as examples.
Furthermore, in that folder contains a notebook for interacting with the [Johnson and Cook (1983)]([johnsonCONSTITUTIVEMODELDATA1983a](@cite)) model.
This [`JohnsonCook.jl`](../../notebooks/JohnsonCook.jl) notebook `include`s a neighboring script ([`JohnsonCook-functions.jl`](../../notebooks/JohnsonCook-functions.jl)) which defines all the necessary types and functions if one wants to create or modify a viscoplasticity model for interaction or calibration as in a Pluto notebook.
The empirical Johnson-Cook model is included as an example given that its widespread recognition might minimize the learning curve for encoding a material model into a Julian implementation appropriate for calibration with `Optimization.jl` or interaction in a Pluto notebook.

## References
```@bibliography
```
