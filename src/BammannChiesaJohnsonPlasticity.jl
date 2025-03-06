module BammannChiesaJohnsonPlasticity


using Reexport
@reexport using ContinuumMechanicsBase


using Tensors # : *, ⊡, sqrt, dev
export δ, vonMises, symmetricmagnitude, symmetricvonMises
δ(i, j) = i == j ? 1.0 : 0.0 # helper function
vonMises(x::SecondOrderTensor) = (s = dev(x); sqrt(3.0/2.0 * s ⊡ s))
symmetricmagnitude(tensor::Vector{<:Real}) = √( sum(tensor[1:3] .^ 2.) + 2sum(tensor[4:6] .^ 2.) )

function symmetricvonMises(tensor::Union{Vector{<:Real}, SubArray{<:Real}}) # ::AbstractFloat
    σvM = sum(map(x->x^2., [tensor[1] - tensor[2], tensor[2] - tensor[3], tensor[3] - tensor[1]])) + (
        6sum(map(x->x^2., [tensor[4], tensor[5], tensor[6]])))
    return √(σvM / 2.)
end

symmetricvonMises(tensor::Matrix{<:Real})::Vector{AbstractFloat} = map(symmetricvonMises, eachcol(tensor))
symmetricvonMises(tensor) = map(symmetricvonMises, eachcol(tensor))


export AbstractBCJ
abstract type AbstractBCJ <: ContinuumMechanicsBase.AbstractMaterialModel end


include("BCJMetal.jl")
include("Bammann1990Modeling.jl")
include("DK.jl")


end
