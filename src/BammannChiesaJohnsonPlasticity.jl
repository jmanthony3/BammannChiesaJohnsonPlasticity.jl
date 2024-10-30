module BammannChiesaJohnsonPlasticity

include("BCJ_solver.jl")
export BCJMetal
export Bammann1990Modeling
export DK
export ISVMetal
export KinematicHardening
export IsotropicHardening
export Damage
export symmetricmagnitude
export symmetricvonMises
export BCJMetalStrainControl
export BCJMetalCurrentConfiguration
export BCJMetalConfigurationHistory
export copyto!
export record!
export bcjmetalreferenceconfiguration
export solve!

include("BCJCalibratinatorJohnsonCookExt.jl")
export JC
export JCStrainControl
export JCCurrentConfiguration
export JCConfigurationHistory
export record!
export jcreferenceconfiguration
export solve!

end
