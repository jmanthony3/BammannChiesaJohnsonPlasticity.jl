module PlasticityCalibratinatorExt

using PlasticityBase
using PlasticityCalibratinator

using BammannChiesaJohnsonPlasticity

include("PlasticityCalibratinatorBCJMetalExt.jl")
include("PlasticityCalibratinatorBammann1990ModelingExt.jl")
include("PlasticityCalibratinatorDKExt.jl")

end