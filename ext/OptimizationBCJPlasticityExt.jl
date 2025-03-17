module OptimizationBCJPlasticityExt

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays
using Optimization, OptimizationOptimJL, LossFunctions

export parameters, parameter_bounds, MaterialOptimizationProblem

include("Metals.jl")

end # end of module
