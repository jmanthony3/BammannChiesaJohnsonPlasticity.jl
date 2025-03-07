module BammannChiesaJohnsonPlasticity


using Reexport
@reexport using ContinuumMechanicsBase


export AbstractBCJ
abstract type AbstractBCJ <: ContinuumMechanicsBase.AbstractMaterialModel end


include("BCJMetal.jl")
include("Bammann1990Modeling.jl")
include("DK.jl")


end
