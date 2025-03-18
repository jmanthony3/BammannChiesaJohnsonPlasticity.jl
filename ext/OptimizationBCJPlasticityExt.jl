module OptimizationBCJPlasticityExt

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays
using Optimization, LossFunctions

export parameters, parameter_bounds, MaterialOptimizationProblem

include("Metals.jl")

end # end of module
