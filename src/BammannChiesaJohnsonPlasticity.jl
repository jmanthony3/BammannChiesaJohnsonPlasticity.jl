module BammannChiesaJohnsonPlasticity



export triu_vec
export δ, vonMises, symmetricmagnitude, symmetricvonMises
export AbstractBCJModel, AbstractBCJTest
export parameters, parameter_bounds, BCJPlasticityProblem


using ContinuumMechanicsBase
using DocStringExtensions
using LinearAlgebra
using Tensors # : *, ⊡, sqrt, dev



# helper functions for working with symmetric tensors
"Get symmetric, second-rank tensor as flat vector."
triu_vec(A::SymmetricTensor{2, 3, <:AbstractFloat}) = A[triu!(trues(size(A)))]

"Kroenecker's Delta"
δ(i, j) = i == j ? 1.0 : 0.0 # helper function

"von Mises (equivalent) scalar for second rank tensor. Internally uses deviatoric."
vonMises(x::SecondOrderTensor) = (s = dev(x); sqrt(3.0/2.0 * s ⊡ s))

"Calculate scalar magnitude for flat vector of symmetric tensor: e. g. A̲̲ ≡⃗ [A₁₁, A₂₂, A₃₃, A₁₂, A₂₃, A₁₃]"
symmetricmagnitude(tensor::Vector{<:Real}) = √( sum(tensor[1:3] .^ 2.) + 2sum(tensor[4:6] .^ 2.) )

"von Mises (equivalent) scalar for flat vector of symmetric tensor: e. g. A̲̲ ≡⃗ [A₁₁, A₂₂, A₃₃, A₁₂, A₂₃, A₁₃]"
function symmetricvonMises(tensor::Union{Vector{<:Real}, SubArray{<:Real}}) # ::AbstractFloat
    σvM = sum(map(x->x^2., [tensor[1] - tensor[2], tensor[2] - tensor[3], tensor[3] - tensor[1]])) + (
        6sum(map(x->x^2., [tensor[4], tensor[5], tensor[6]])))
    return √(σvM / 2.)
end
# symmetricvonMises(tensor::Matrix{<:Real})::Vector{AbstractFloat} = map(symmetricvonMises, eachcol(tensor))
symmetricvonMises(tensor) = map(symmetricvonMises, eachcol(tensor))



# module specific codes from parent, [CMB.jl](https://github.com/TRACER-LULab/ContinuumMechanicsBase.jl.git)
"Define parent type for all BCJ-variant models."
abstract type AbstractBCJModel      <: ContinuumMechanicsBase.AbstractMaterialModel end

"Define parent type for all BCJ-variant tests."
abstract type AbstractBCJTest       <: ContinuumMechanicsBase.AbstractMaterialTest end


include("Metals.jl") # `include`s for metal-specific models



end # that's all folks!
