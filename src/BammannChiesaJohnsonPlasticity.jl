module BammannChiesaJohnsonPlasticity


using Reexport
@reexport using ContinuumMechanicsBase


using Tensors # : *, ⊡, sqrt, dev
export δ, vonMises
δ(i, j) = i == j ? 1.0 : 0.0 # helper function
vonMises(x) = (s = dev(x); sqrt(3.0/2.0 * s ⊡ s))


export AbstractBCJ
abstract type AbstractBCJ <: ContinuumMechanicsBase.AbstractMaterialModel end


include("BCJMetal.jl")
include("Bammann1990Modeling.jl")
include("DK.jl")


end
