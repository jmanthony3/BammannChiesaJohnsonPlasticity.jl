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


include("BCJMetals.jl") # `include`s for metal-specific models


# function parameters(::M) where {M<:ContinuumMechanicsBase.AbstractMaterialModel} end

# function parameter_bounds(::M, ::Any) where {M<:ContinuumMechanicsBase.AbstractMaterialModel}
#     lb = nothing
#     ub = nothing
#     return (lb = lb, ub = ub)
# end

# function parameter_bounds(
#             ψ       ::M,
#             tests   ::Vector{Any},
#         ) where {M<:ContinuumMechanicsBase.AbstractMaterialModel}
#     bounds = map(Base.Fix1(parameter_bounds, ψ), tests)
#     lbs = getfield.(bounds, :lb)
#     ubs = getfield.(bounds, :ub)
#     if !(eltype(lbs) <: Nothing)
#         lb_ps = fieldnames(eltype(lbs))
#         lb = map(p -> p .=> maximum(getfield.(lbs, p)), lb_ps) |> NamedTuple
#     else
#         lb = nothing
#     end
#     if !(eltype(ubs) <: Nothing)
#         ub_ps = fieldnames(eltype(ubs))
#         ub = map(p -> p .=> minimum(getfield.(ubs, p)), ub_ps) |> NamedTuple
#     else
#         ub = nothing
#     end
#     return (lb = lb, ub = ub)
# end

# """
# $(SIGNATURES)

# Creates an `OptimizationProblem` for use in [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) to find the optimal parameters.

# # Arguments:
# - `ψ`: material model to use
# - `test` or `tests`: A single or vector of hyperelastics tests to use when fitting the parameters
# - `u₀`: Initial guess for parameters
# - `ps`: Any additional parameters for calling predict
# - `adb`: Select differentiation type from [`ADTypes.jl`](https://github.com/SciML/ADTypes.jl). The type is automatically applied to the type of AD applied to the Optimization Problem also.
# - `loss`: Loss function from [`LossFunctions.jl`](https://github.com/JuliaML/LossFunctions.jl)
# """
# function BCJPlasticityProblem end



end # that's all folks!
