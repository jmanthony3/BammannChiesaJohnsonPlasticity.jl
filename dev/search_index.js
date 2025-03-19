var documenterSearchIndex = {"docs":
[{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/#OptimizationBCJPlasticityExt","page":"Optimization.jl","title":"OptimizationBCJPlasticityExt","text":"","category":"section"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"CurrentModule = Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)","category":"page"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"This extension is meant to provide additional functionality to the base package with the intent of performing optimization studies with the Optimization.jl and LossFunctions.jl packages to calibrate material model constants against experimental data. This interface with Optimization.jl is largely accomplished by overloading three functions from ContinuumMechanicsBase.jl with the newly defined types in the BCJPlasticity.jl package:","category":"page"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"parameters,\nparameter_bounds, and\nMaterialOptimizationProblem","category":"page"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"Most significant of these overloads is that for parameter_bounds wherein the default value for lb of the returned tuple, which represents the lower bounds of possible values during optimization, is constrained to be 0 infty) to maintain physical admissibility and self-consistency in the constitutive and ISV rate equations.","category":"page"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/#Functions","page":"Optimization.jl","title":"Functions","text":"","category":"section"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"Two functions are overloaded for any sub-type of AbstractBCJModel and AbstractBCJTest: parameter_bounds and MaterialOptimizationProblem. The overload for parameter_bounds, by default, calls on parameters which must be overloaded for the appropriate BCJ model.","category":"page"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"Modules = [Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)]\nOrder   = [:function]\nPages   = [\"OptimizationBCJPlasticityExt.jl\"]","category":"page"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/#ContinuumMechanicsBase.MaterialOptimizationProblem-Tuple{AbstractBCJModel, AbstractBCJTest, Vararg{Any, 4}}","page":"Optimization.jl","title":"ContinuumMechanicsBase.MaterialOptimizationProblem","text":"Dispatch for BCJ-specific types and functions.\n\n\n\n\n\n","category":"method"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/#ContinuumMechanicsBase.parameter_bounds-Tuple{AbstractBCJModel, AbstractBCJTest}","page":"Optimization.jl","title":"ContinuumMechanicsBase.parameter_bounds","text":"Set lower bounds to zero for physical admissibility.\n\n\n\n\n\n","category":"method"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/#Metals","page":"Optimization.jl","title":"Metals","text":"","category":"section"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"Modules = [Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)]\nOrder   = [:function]\nPages   = [\"Metals.jl\",\n    \"Bammann1990Modeling.jl\"]","category":"page"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/#ContinuumMechanicsBase.parameters-Tuple{Bammann1990Modeling}","page":"Optimization.jl","title":"ContinuumMechanicsBase.parameters","text":"Constants for temperature equations from Bammann (1990). Note: though not explicitly listed in paper, temperature equations h = C₁₅ * exp(C₁₆ * θ) and H = C₁₇ * exp(C₁₈ * θ) are included (c. f. DYNA3D User Manual (1993)).\n\n\n\n\n\n","category":"method"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/#Index","page":"Optimization.jl","title":"Index","text":"","category":"section"},{"location":"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt/","page":"Optimization.jl","title":"Optimization.jl","text":"Modules = [Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)]\nOrder   = [:type, :function]\nPages   = [@__FILE__]","category":"page"},{"location":"base/Bammann1990Modeling/#Bammann1990Modeling","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"","category":"section"},{"location":"base/Bammann1990Modeling/","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"CurrentModule = BammannChiesaJohnsonPlasticity","category":"page"},{"location":"base/Bammann1990Modeling/","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"The following equations are those employed in the Bammann (1990) paper that are implemented in the Bammann1990Modeling type, constructor, and that method of the kernel function, update associate with this type. As such, the internal equations use the same nomenclature for plastic strain rate, mathbfD^p equiv dotepsilon^(p); the second-rank, deviatoric tensors for Cauchy stress, undersetsimsigma and kinematic hardening, undersetsimalpha which is an ISV; the scalar isotropic hardening, kappa, and the constants associated with the dynamic and static recovery temperature equations. Of major importance is that, although not explicitly listed in the publication, the equations for h and H are included in this implementation (c. f. DYNA3D User Manual (1993)).","category":"page"},{"location":"base/Bammann1990Modeling/","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"beginaligned\n     plastic strain rate\n    mathbfD^p = f(theta)sinhleft frac undersetsimxi - kappa - Y(theta)  V(theta)  rightfracundersetsimxiundersetsimxitext let undersetsimxi = undersetsimsigma - undersetsimalpha \n     kinematic hardening\n    dotundersetsimalpha = hmu(theta)mathbfD^p - r_d(theta)mathbfD^p + r_s(theta)undersetsimalphaundersetsimalpha \n     isotropic hardening\n    dotkappa = Hmu(theta)mathbfD^p - R_d(theta)mathbfD^p + R_s(theta)kappa^2 \n     flow rule\n    F = undersetsimsigma - undersetsimalpha - kappa - beta(mathbfD^p theta) \n     initial yield stress beta\n    beta(mathbfD^p theta) = Y(theta) + V(theta)sinh^-1left(fracmathbfD^pf(theta)right) \n    V(theta)       = C_ 1 expleft( -frac C_ 2  theta  right) \n    Y(theta)       = C_ 3 expleft(  frac C_ 4  theta  right) \n    f(theta)       = C_ 5 expleft( -frac C_ 6  theta  right) \n    r_d(theta)   = C_ 7 expleft( -frac C_ 8  theta  right) \n    r_s(theta)   = C_ 9 expleft( -frac C_10  theta  right) \n    R_d(theta)   = C_11 expleft( -frac C_12  theta  right) \n    R_s(theta)   = C_13 expleft( -frac C_14  theta  right) \n    h(theta)       = C_15 expleft(  frac C_16  theta  right) \n    H(theta)       = C_17 expleft(  frac C_18  theta  right) \nendaligned","category":"page"},{"location":"base/Bammann1990Modeling/#Types","page":"Bammann1990Modeling","title":"Types","text":"","category":"section"},{"location":"base/Bammann1990Modeling/","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"Modules = [BammannChiesaJohnsonPlasticity, ContinuumMechanicsBase]\nOrder   = [:type]\nPages   = [\"Bammann1990Modeling.jl\"]","category":"page"},{"location":"base/Bammann1990Modeling/#BammannChiesaJohnsonPlasticity.Bammann1990Modeling","page":"Bammann1990Modeling","title":"BammannChiesaJohnsonPlasticity.Bammann1990Modeling","text":"Structure for viscoplasticity model with loading conditions and material properties. Here, uses the effective strain rate based on applied strain rate and loading direction.\n\n\n\n\n\n","category":"type"},{"location":"base/Bammann1990Modeling/#BammannChiesaJohnsonPlasticity.Bammann1990Modeling-Tuple{BCJMetalStrainControl, AbstractFloat}","page":"Bammann1990Modeling","title":"BammannChiesaJohnsonPlasticity.Bammann1990Modeling","text":"Bammann1990Modeling(Ω, μ)\n\n\nOuter constructor for loading conditions and material properties which assumes a Poisson's ratio of 0.5. Here, μ is the shear modulus.\n\n\n\n\n\n","category":"method"},{"location":"base/Bammann1990Modeling/#Functions","page":"Bammann1990Modeling","title":"Functions","text":"","category":"section"},{"location":"base/Bammann1990Modeling/","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"Modules = [BammannChiesaJohnsonPlasticity, ContinuumMechanicsBase]\nOrder   = [:function]\nPages   = [\"Bammann1990Modeling.jl\"]","category":"page"},{"location":"base/Bammann1990Modeling/#BammannChiesaJohnsonPlasticity.update-Tuple{Bammann1990Modeling, Vararg{Any, 6}}","page":"Bammann1990Modeling","title":"BammannChiesaJohnsonPlasticity.update","text":"Using the equations and constants from Bammann (1990), this kernel function maps the current material state and ISVs onto the next configuration. Note: though not explicitly listed in paper, temperature equations h = C₁₅ * exp(C₁₆ * θ) and H = C₁₇ * exp(C₁₈ * θ) are included (c. f. DYNA3D User Manual (1993)).\n\n\n\n\n\n","category":"method"},{"location":"base/Bammann1990Modeling/#References","page":"Bammann1990Modeling","title":"References","text":"","category":"section"},{"location":"base/Bammann1990Modeling/","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"D. J. Bammann. Modeling Temperature and Strain Rate Dependent Large Deformations of Metals. Applied Mechanics Reviews 43, S312-S319 (1990).\n\n\n\nR. G. Whirley and B. Engelmann. DYNA3D: A Nonlinear, Explicit, Three-Dimensional Finite Element Code for Solid and Structural Mechanics, User Manual. Revision 1 (Lawrence Livermore National Lab.(LLNL), Livermore, CA (United States), 1993).\n\n\n\n","category":"page"},{"location":"base/Bammann1990Modeling/#Index","page":"Bammann1990Modeling","title":"Index","text":"","category":"section"},{"location":"base/Bammann1990Modeling/","page":"Bammann1990Modeling","title":"Bammann1990Modeling","text":"Modules = [BammannChiesaJohnsonPlasticity, ContinuumMechanicsBase]\nOrder   = [:type, :function]\nPages   = [@__FILE__]","category":"page"},{"location":"primer/#Primer","page":"Primer","title":"Primer","text":"","category":"section"},{"location":"primer/","page":"Primer","title":"Primer","text":"Bammann-Chiesa-Johnson (BCJ) plasticity models the inelastic deformation of solid materials at a given temperature and strain rate with Internal State Variables (ISVs). The BCJ-plasticity model satisfies the Clausius-Duhem inequality by placing thermodynamic constraints on ISVs and their rate equations (Coleman and Gurtin (1967)) which are motivated by equations for dislocation mechanics. The model assumes that for any continuum mechanics configuration each material point at each moment in time satisfies some equilibrium criterion for a triplet point composed of the deformation gradient, temperature, and current set of ISVs. The ISVs in the current configuration intrinsically represent their evolution because of the material history afforded by the rate equations. These ISVs and their rate equations are used to increment the current configuration from the previous configuration by a trial guess that the material deforms elastically and correcting this trial guess for any plastic deformation accrued in the increment: this is the scheme of the Radial-Return Method from Krieg and Krieg (1977). Using ISVs to model inelastic deformation is not unique to BCJ-plasticity model. For a broader review of modeling inelastic deformation with ISVs, as well as, a history of the BCJ-plasticity model, the reader is directed Dr. Horstemeyer's 2010 review paper [4].","category":"page"},{"location":"primer/#References","page":"Primer","title":"References","text":"","category":"section"},{"location":"primer/","page":"Primer","title":"Primer","text":"B. D. Coleman and M. E. Gurtin. Thermodynamics with Internal State Variables. The Journal of Chemical Physics 47, 597–613 (1967).\n\n\n\nR. D. Krieg and D. B. Krieg. Accuracies of Numerical Solution Methods for the Elastic-Perfectly Plastic Model. Journal of Pressure Vessel Technology 99, 510–515 (1977).\n\n\n\nM. F. Horstemeyer and D. J. Bammann. Historical Review of Internal State Variable Theory for Inelasticity. International Journal of Plasticity 26, 1310–1334 (2010).\n\n\n\n","category":"page"},{"location":"base/BammannChiesaJohnsonPlasticity/#Base","page":"Base Package","title":"Base","text":"","category":"section"},{"location":"base/BammannChiesaJohnsonPlasticity/","page":"Base Package","title":"Base Package","text":"CurrentModule = BammannChiesaJohnsonPlasticity","category":"page"},{"location":"base/BammannChiesaJohnsonPlasticity/","page":"Base Package","title":"Base Package","text":"This package adds overloads and sub-types from the ContinuumMechanicsBase.jl package for the BCJ-plasticity model. Most importantly, the AbstractBCJModel and AbstractBCJTest sub-types are defined first such that subsequent types and functions might dispatch under these two types. Some helper functions are also defined for calculating the mathematical norm from a vector representing the upper triangle of a symmetric tensor and the von Mises equivalent of a symmetric tensor. The elements of this vector for the upper triangle of a symmetric vector are in row-major order. That is, the upper triangle of a symmetric, second-rank tensor may be represented as a vector: e. g. mathbfA vecequiv A_11 A_12 A_13 A_22 A_23 A_33. Currently, the package includes support for the following categories of materials:","category":"page"},{"location":"base/BammannChiesaJohnsonPlasticity/","page":"Base Package","title":"Base Package","text":"Modules = [BammannChiesaJohnsonPlasticity]\nPages   = [\n    \"base/Metals.md\"]\nDepth   = 1","category":"page"},{"location":"base/BammannChiesaJohnsonPlasticity/#Types","page":"Base Package","title":"Types","text":"","category":"section"},{"location":"base/BammannChiesaJohnsonPlasticity/","page":"Base Package","title":"Base Package","text":"Modules = [BammannChiesaJohnsonPlasticity]\nOrder   = [:type]\nPages   = [\"BammannChiesaJohnsonPlasticity.jl\"]","category":"page"},{"location":"base/BammannChiesaJohnsonPlasticity/#BammannChiesaJohnsonPlasticity.AbstractBCJModel","page":"Base Package","title":"BammannChiesaJohnsonPlasticity.AbstractBCJModel","text":"Parent type for all BCJ-variant models.\n\n\n\n\n\n","category":"type"},{"location":"base/BammannChiesaJohnsonPlasticity/#BammannChiesaJohnsonPlasticity.AbstractBCJTest","page":"Base Package","title":"BammannChiesaJohnsonPlasticity.AbstractBCJTest","text":"Parent type for all BCJ-variant tests.\n\n\n\n\n\n","category":"type"},{"location":"base/BammannChiesaJohnsonPlasticity/#Functions","page":"Base Package","title":"Functions","text":"","category":"section"},{"location":"base/BammannChiesaJohnsonPlasticity/","page":"Base Package","title":"Base Package","text":"Modules = [BammannChiesaJohnsonPlasticity]\nOrder   = [:function]\nPages   = [\"BammannChiesaJohnsonPlasticity.jl\"]","category":"page"},{"location":"base/BammannChiesaJohnsonPlasticity/#BammannChiesaJohnsonPlasticity.norm_symvec-Tuple{Vector{<:Real}}","page":"Base Package","title":"BammannChiesaJohnsonPlasticity.norm_symvec","text":"Calculate scalar magnitude for flat vector of symmetric tensor: e. g. A̲̲ ≡⃗ [A₁₁, A₁₂, A₁₃, A₂₂, A₂₃, A₃₃]\n\n\n\n\n\n","category":"method"},{"location":"base/BammannChiesaJohnsonPlasticity/#BammannChiesaJohnsonPlasticity.update-Tuple{AbstractBCJModel, Vararg{Any}}","page":"Base Package","title":"BammannChiesaJohnsonPlasticity.update","text":"update(ψ, args; kwargs...)\n\n\nGiven viscoplasticity model and the current material state, update to the next material state.\n\n\n\n\n\n","category":"method"},{"location":"base/BammannChiesaJohnsonPlasticity/#BammannChiesaJohnsonPlasticity.vonMises-Tuple{Any}","page":"Base Package","title":"BammannChiesaJohnsonPlasticity.vonMises","text":"von Mises (equivalent) scalar for symmetric tensor.\n\n\n\n\n\n","category":"method"},{"location":"base/BammannChiesaJohnsonPlasticity/#BammannChiesaJohnsonPlasticity.vonMises-Tuple{Union{SubArray{<:Real}, Vector{<:Real}}}","page":"Base Package","title":"BammannChiesaJohnsonPlasticity.vonMises","text":"von Mises (equivalent) scalar for flat vector of symmetric tensor: e. g. A̲̲ ≡⃗ [A₁₁, A₁₂, A₁₃, A₂₂, A₂₃, A₃₃]\n\n\n\n\n\n","category":"method"},{"location":"base/BammannChiesaJohnsonPlasticity/#Index","page":"Base Package","title":"Index","text":"","category":"section"},{"location":"base/BammannChiesaJohnsonPlasticity/","page":"Base Package","title":"Base Package","text":"Modules = [BammannChiesaJohnsonPlasticity]\nOrder   = [:type, :function]\nPages   = [@__FILE__]","category":"page"},{"location":"base/Metals/#Metals","page":"Metals","title":"Metals","text":"","category":"section"},{"location":"base/Metals/","page":"Metals","title":"Metals","text":"CurrentModule = BammannChiesaJohnsonPlasticity","category":"page"},{"location":"base/Metals/","page":"Metals","title":"Metals","text":"Modules = [BammannChiesaJohnsonPlasticity]\nPages   = [\n    \"base/Bammann1990Modeling.md\"]\nDepth   = 1","category":"page"},{"location":"base/Metals/#Types","page":"Metals","title":"Types","text":"","category":"section"},{"location":"base/Metals/","page":"Metals","title":"Metals","text":"Modules = [BammannChiesaJohnsonPlasticity]\nOrder   = [:type]\nPages   = [\"Metals.jl\"]","category":"page"},{"location":"base/Metals/#BammannChiesaJohnsonPlasticity.AbstractBCJMetalModel","page":"Metals","title":"BammannChiesaJohnsonPlasticity.AbstractBCJMetalModel","text":"Sub-type of AbstractBCJModel for BCJ-models specific to metals.\n\n\n\n\n\n","category":"type"},{"location":"base/Metals/#BammannChiesaJohnsonPlasticity.AbstractBCJMetalTest","page":"Metals","title":"BammannChiesaJohnsonPlasticity.AbstractBCJMetalTest","text":"Sub-type of AbstractBCJTest for BCJ-models specific to metals.\n\n\n\n\n\n","category":"type"},{"location":"base/Metals/#BammannChiesaJohnsonPlasticity.BCJMetalDataEntry","page":"Metals","title":"BammannChiesaJohnsonPlasticity.BCJMetalDataEntry","text":"Store vectors of strain, ϵ, and stress, σ data.\n\n\n\n\n\n","category":"type"},{"location":"base/Metals/#BammannChiesaJohnsonPlasticity.BCJMetalStrainControl","page":"Metals","title":"BammannChiesaJohnsonPlasticity.BCJMetalStrainControl","text":"Stucture for strain-controlled loadings of metals for temperature, θ; strain rate, ϵ̇; final strain, ϵₙ; number of loading increments, N; loading direction, loaddir ∈ {:tension, :compression, :torsion}.\n\n\n\n\n\n","category":"type"},{"location":"base/Metals/#Functions","page":"Metals","title":"Functions","text":"","category":"section"},{"location":"base/Metals/","page":"Metals","title":"Metals","text":"Modules = [BammannChiesaJohnsonPlasticity]\nOrder   = [:function]\nPages   = [\"Metals.jl\"]","category":"page"},{"location":"base/Metals/#Index","page":"Metals","title":"Index","text":"","category":"section"},{"location":"base/Metals/","page":"Metals","title":"Metals","text":"Modules = [BammannChiesaJohnsonPlasticity]\nOrder   = [:type, :function]\nPages   = [@__FILE__]","category":"page"},{"location":"#BammannChiesaJohnsonPlasticity","page":"Home","title":"BammannChiesaJohnsonPlasticity","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BammannChiesaJohnsonPlasticity.","category":"page"},{"location":"#Purpose","page":"Home","title":"Purpose","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The BammannChiesaJohnsonPlasticity.jl package is a Julian implementation of Bammann-Chiesa-Johnson (BCJ) plasticity from select publications. The goal of this package is to provide a repository and interface for the BCJ-plasticity model and its equations to implement and interact with those equations from publications, as well as, give opportunity for developing or improving the model. The BCJ model has been applied to many materials over the years; however, only those models from a select few papers for modeling the inelastic deformation of metals are currently defined in this package. It is a hope of the authors that published models expanding, refining, or applying the BCJ model onto a wider selection of materials might be included in this package along with the model constants employed in those publications for the appropriate chemical system. Having such a repository for versions of the BCJ-plasticity model, variations of coupling BCJ-plasticity with other ISV models (damage, recrystallization, etcetera), and the relevant constants maintains transparency and reproducibility of publications. Many publications discuss the physics and mechanics of the BCJ-plasticity model, but a small primer is included in the next section for those unfamiliar with the model. At this time, all other sections after the primer will only include comments on the behavior of the package and not the limitations of the model as this is discussed in the literature.","category":"page"},{"location":"#Methodology","page":"Home","title":"Methodology","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This BCJPlasticity.jl package employs and overloads the ContinuumMechanicsBase.jl package for the types and functions common to continuum mechanics modeling. Currently, the package acts as a point simulator to update the material states and assumes a Poisson's ratio of 0.5. Although the package can be used out-of-the-box to plot the specified BCJ model with the appropriate loading conditions, material properties, and model constants, the package can be extended to calibrate the constants against experimental data via the Optimization.jl package. For an example of interacting with and calibrating BCJ model constants, the reader is referred to the notebooks folder for a set of Pluto.jl notebooks as examples. Furthermore, in that folder contains a notebook for interacting with the Johnson and Cook (1983). This JohnsonCook.jl notebook includes a neighboring script (JohnsonCook-functions.jl) which defines all the necessary types and functions if one wants to create or modify a viscoplasticity model for interaction or calibration as in a Pluto notebook. The empirical Johnson-Cook model is included as an example given that its widespread recognition might minimize the learning curve for encoding a material model into a Julian implementation appropriate for calibration with Optimization.jl or interaction in a Pluto notebook.","category":"page"},{"location":"#Nomenclature","page":"Home","title":"Nomenclature","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As much as possible, equations listed in the documentation for any version of the BCJ-plasticity model will attempt to match the form of the equations from the cited publication. For consistency, some cosmetics are applied in the documentation to maintain clarity as to the rank of tensors. For English letters, first and second rank tensors will be denoted as lower and uppercase, bold-face symbols, respectively: e. g. mathbfx is a first rank tensor while mathbfA is second rank. For Greek letters, Einsteinian notation, instead of bold type facing, is used to connote tensors: e. g. undersetsimsigma could be the second rank tensor for stress. Many functions used throughout this package, as well as its overloads for ContinuumMechanicsBase.jl, assign the prescribed BCJ-plasticity model to the variable ψ and the applied boundary conditions assigned to Ω.","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"G. R. Johnson and W. H. Cook. A CONSTITUTIVE MODEL AND DATA FOR METALS SUBJECTED TO LARGE STRAINS, HIGH STRAIN RATES AND HIGH TEMPERATURES. In: Proceedings of the 7th International Symposium on Ballistics (The Hague, Netherlands, 1983); pp. 541–547.\n\n\n\nB. D. Coleman and M. E. Gurtin. Thermodynamics with Internal State Variables. The Journal of Chemical Physics 47, 597–613 (1967).\n\n\n\nR. D. Krieg and D. B. Krieg. Accuracies of Numerical Solution Methods for the Elastic-Perfectly Plastic Model. Journal of Pressure Vessel Technology 99, 510–515 (1977).\n\n\n\nM. F. Horstemeyer and D. J. Bammann. Historical Review of Internal State Variable Theory for Inelasticity. International Journal of Plasticity 26, 1310–1334 (2010).\n\n\n\nD. J. Bammann. Modeling Temperature and Strain Rate Dependent Large Deformations of Metals. Applied Mechanics Reviews 43, S312-S319 (1990).\n\n\n\nR. G. Whirley and B. Engelmann. DYNA3D: A Nonlinear, Explicit, Three-Dimensional Finite Element Code for Solid and Structural Mechanics, User Manual. Revision 1 (Lawrence Livermore National Lab.(LLNL), Livermore, CA (United States), 1993).\n\n\n\n","category":"page"},{"location":"extensions/#Extensions","page":"Extending Functionality","title":"Extensions","text":"","category":"section"},{"location":"extensions/","page":"Extending Functionality","title":"Extending Functionality","text":"The functionality of the base package can be extended upon loading the other packages. Currently, the following extensions can be loaded for additional interaction with the BCJ-plasticity model, and the reader is directed to those pages for the extension of interest. Future work may include an extension for the Ferrite.jl package to perform finite element analyses in pure Julia with the BCJ-plasticity model.","category":"page"},{"location":"extensions/","page":"Extending Functionality","title":"Extending Functionality","text":"Modules = [\n    Base.get_extension(BammannChiesaJohnsonPlasticity, :OptimizationBCJPlasticityExt)]\nPages   = [\n    \"ext/OptimizationBCJPlasticityExt.jl/OptimizationBCJPlasticityExt.md\"]\nDepth   = 1","category":"page"}]
}
