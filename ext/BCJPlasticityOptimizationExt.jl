module BCJPlasticityOptimizationExt

using BammannChiesaJohnsonPlasticity

using ComponentArrays
using Optimization, OptimizationOptimJL, LossFunctions

export parameters, parameter_bounds, BCJPlasticityProblem

include("BCJMetals.jl")

end # end of module
