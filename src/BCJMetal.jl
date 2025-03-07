using ContinuumMechanicsBase

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



export BCJMetal, BCJMetalStrainControl, AbstractBCJMetalTest, BCJMetalDataEntry, BCJMetalUniaxialTest
export parameter_bounds

abstract type BCJMetal                      <: AbstractBCJ end
# abstract type ISVMetal{T<:BCJMetal} end
# abstract type ISVMetalKinematicHardening{T} <: ISVMetal{T} end # α__
# abstract type ISVMetalIsotropicHardening{T} <: ISVMetal{T} end # κ
# abstract type ISVMetalDamage{T}             <: ISVMetal{T} end # ϕ

struct BCJMetalStrainControl{T1<:Integer, T2<:AbstractFloat} <: ContinuumMechanicsBase.AbstractMaterialTest
    θ           ::T2                # applied temperature
    ϵ_dot       ::T2                # applied strain rate
    ϵₙ          ::T2                # final strain
    N           ::T1                # number of strain increments
    loadtype    ::Symbol            # load type (:tension, :compression, :torsion)
    # constants   ::OrderedDict{String, T2}  # material constants
end

abstract type AbstractBCJMetalTest{T, S} <: ContinuumMechanicsBase.AbstractMaterialTest end

using ComponentArrays, StructArrays
struct BCJMetalDataEntry{T, S}
    λ::Vector{T}
    s::Vector{S}
end

struct BCJMetalUniaxialTest{T, S} <: AbstractBCJMetalTest{T, S}
    data::StructVector
    name::String
    """
    $(SIGNATURES)

    Creates an object storing results from a uniaxial test of a hyperelatic  material.

    # Arguments:
    - `λ₁`: Vector of uniaxial stretches
    - `s₁`: Vector of experiemntal stresses (optional)
    - `name`: string for the name of the test
    - `incompressible`: `true` if the material can be assumed to be incompressible.
    """
    function BCJMetalUniaxialTest(λ₁, s₁; name, incompressible = true)
        @assert length(λ₁) == length(s₁) "Inputs must be the same length"
        if incompressible
            # λ₂ = λ₃ = @. sqrt(1 / λ₁)
            λ₂ = λ₃ = @. -0.499λ₁
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁))
        data = StructArray{BCJMetalDataEntry}((λ, s))
        new{eltype(eltype(λ)),eltype(eltype(s))}(data, name)
    end
    function BCJMetalUniaxialTest(λ₁; name, incompressible = true)
        if incompressible
            # λ₂ = λ₃ = @. sqrt(1 / λ₁)
            λ₂ = λ₃ = @. -0.499λ₁
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(Vector{eltype(λ₁)}(undef, length(λ₁))))
        data = StructArray{BCJMetalDataEntry}((λ, s))
        new{eltype(eltype(λ)),eltype(eltype(s))}(data, name)
    end
end

function parameter_bounds(::ContinuumMechanicsBase.AbstractMaterialModel, ::Any)
    lb = nothing
    ub = nothing
    return (lb = lb, ub = ub)
end

function parameter_bounds(
            ψ       ::ContinuumMechanicsBase.AbstractMaterialModel,
            tests   ::Vector{Any},
        )
    bounds = map(Base.Fix1(parameter_bounds, ψ), tests)
    lbs = getfield.(bounds, :lb)
    ubs = getfield.(bounds, :ub)
    if !(eltype(lbs) <: Nothing)
        lb_ps = fieldnames(eltype(lbs))
        lb = map(p -> p .=> maximum(getfield.(lbs, p)), lb_ps) |> NamedTuple
    else
        lb = nothing
    end
    if !(eltype(ubs) <: Nothing)
        ub_ps = fieldnames(eltype(ubs))
        ub = map(p -> p .=> minimum(getfield.(ubs, p)), ub_ps) |> NamedTuple
    else
        ub = nothing
    end
    return (lb = lb, ub = ub)
end