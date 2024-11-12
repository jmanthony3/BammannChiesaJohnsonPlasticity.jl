module BammannChiesaJohnsonPlasticity

using PlasticityBase

export kernel
export nextloadingphase

abstract type BCJ <: AbstractPlasticity end

export BCJ

include("BCJMetal.jl")
export BCJMetal
export ISVMetal
export ISVMetalKinematicHardening
export ISVMetalIsotropicHardening
export ISVMetalDamage
export BCJMetalStrainControl
export BCJMetalConfigurationCurrent
export BCJMetalConfigurationHistory
export BCJMetalConfigurationTuple
export +
export copyto!
export record!
export symmetricmagnitude
export symmetricvonMises

include("Bammann1990Modeling.jl")
export Bammann1990Modeling
export referenceconfiguration
export solve!

include("DK.jl")
export DK
export referenceconfiguration
export solve!

end
