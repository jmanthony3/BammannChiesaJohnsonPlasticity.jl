module BammannChiesaJohnsonPlasticity



# export triu_vec, δ, vonMises # uncomment when we can work with Tensors.jl (#5)
export norm_symvec, vonMises
export AbstractBCJModel, AbstractBCJTest, update


using ContinuumMechanicsBase
using DocStringExtensions
using LinearAlgebra
# using Tensors # : *, ⊡, sqrt, dev



# # helper functions for working with symmetric tensors
# "Get symmetric, second-rank tensor as flat vector."
# triu_vec(A::SymmetricTensor{2, 3, <:AbstractFloat}) = A[triu!(trues(size(A)))]

# "Kroenecker's Delta"
# δ(i, j) = i == j ? 1.0 : 0.0 # helper function

# "von Mises (equivalent) scalar for second rank tensor. Internally uses deviatoric."
# vonMises(x::SecondOrderTensor) = (s = dev(x); sqrt(3.0/2.0 * s ⊡ s))

"Calculate scalar magnitude for flat vector of symmetric tensor: e. g. A̲̲ ≡⃗ [A₁₁, A₁₂, A₁₃, A₂₂, A₂₃, A₃₃]"
norm_symvec(tensor::Vector{<:Real}) = √( sum(tensor[[1, 4, 6]] .^ 2.0) + 2sum(tensor[[2, 3, 5]] .^ 2.0) )

"von Mises (equivalent) scalar for flat vector of symmetric tensor: e. g. A̲̲ ≡⃗ [A₁₁, A₁₂, A₁₃, A₂₂, A₂₃, A₃₃]"
function vonMises(tensor::Union{Vector{<:Real}, SubArray{<:Real}}) # ::AbstractFloat
    σvM = sum(map(x->x^2.0, [tensor[1] - tensor[4], tensor[4] - tensor[6], tensor[6] - tensor[1]])) + (
        6sum(map(x->x^2.0, [tensor[2], tensor[3], tensor[5]])))
    return √(σvM / 2.)
end
# symmetricvonMises(tensor::Matrix{<:Real})::Vector{AbstractFloat} = map(symmetricvonMises, eachcol(tensor))
"von Mises (equivalent) scalar for symmetric tensor."
vonMises(tensor) = map(vonMises, eachcol(tensor))



# module specific codes from parent, [CMB.jl](https://github.com/TRACER-LULab/ContinuumMechanicsBase.jl.git)
"Parent type for all BCJ-variant models."
abstract type AbstractBCJModel      <: ContinuumMechanicsBase.AbstractMaterialModel end

"Parent type for all BCJ-variant tests."
abstract type AbstractBCJTest       <: ContinuumMechanicsBase.AbstractMaterialTest end

"""
    $(SIGNATURES)

Given viscoplasticity model and the current material state, update to the next material state.
"""
function update(ψ::AbstractBCJModel, args...; kwargs...) end


include("Metals.jl") # `include`s for metal-specific models



end # that's all folks!
