module OptimizationBCJPlasticityExt

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays
using Optimization, OptimizationOptimJL, LossFunctions

export parameters, parameter_bounds, MaterialOptimizationProblem

include("BCJMetals.jl")

end # end of module
