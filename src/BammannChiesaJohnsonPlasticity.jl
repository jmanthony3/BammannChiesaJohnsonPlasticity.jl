module BammannChiesaJohnsonPlasticity

using PlasticityBase

abstract type BCJ <: Plasticity end

export BCJ

include("BCJMetal.jl")
export BCJMetal
export ISVMetal
export KinematicHardening
export IsotropicHardening
export Damage
export symmetricmagnitude
export symmetricvonMises
export BCJMetalStrainControl
export BCJMetalCurrentConfiguration
export BCJMetalConfigurationHistory
export record!

include("Bammann1990Modeling.jl")
export Bammann1990Modeling
export referenceconfiguration
export solve!

include("DK.jl")
export DK
export referenceconfiguration
export solve!

end
