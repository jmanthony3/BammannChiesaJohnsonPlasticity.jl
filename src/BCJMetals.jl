# module Metals



export AbstractBCJMetalModel
export BCJMetalStrainControl
export BCJMetalDataEntry
export BCJMetalUniaxialTest

using ContinuumMechanicsBase
using ComponentArrays, StructArrays



abstract type AbstractBCJMetalModel         <: AbstractBCJModel end
# abstract type ISVMetal{T<:BCJMetal} end
# abstract type ISVMetalKinematicHardening{T} <: ISVMetal{T} end # α__
# abstract type ISVMetalIsotropicHardening{T} <: ISVMetal{T} end # κ
# abstract type ISVMetalDamage{T}             <: ISVMetal{T} end # ϕ
abstract type AbstractBCJMetalTest{T}       <: AbstractBCJTest end

struct BCJMetalStrainControl{T1<:Integer, T2<:AbstractFloat} <: AbstractBCJTest
    θ           ::T2                # applied temperature
    ϵ_dot       ::T2                # applied strain rate
    ϵₙ          ::T2                # final strain
    N           ::T1                # number of strain increments
    loadtype    ::Symbol            # load type (:tension, :compression, :torsion)
end

struct BCJMetalDataEntry{T, S}
    ϵ::Vector{T}
    σ::Vector{S}
end

struct BCJMetalUniaxialTest{T, S} <: AbstractBCJMetalTest{T}
    data::StructVector
    name::String
    """
    $(SIGNATURES)

    Creates an object storing results from a uniaxial test of a hyperelatic  material.

    # Arguments:
    - `ϵ₁`: Vector of experimental, uniaxial strains
    - `σ₁`: Vector of experimental, uniaxial stresses (optional)
    - `name`: string for the name of the test
    - `incompressible`: `true` if the material can be assumed to be incompressible.
    """
    function BCJMetalUniaxialTest(ϵ₁, σ₁; name, incompressible = true)
        @assert length(ϵ₁) == length(σ₁) "Inputs must be the same length"
        if incompressible
            # ϵ₂ = ϵ₃ = @. sqrt(1 / ϵ₁)
            ϵ₂ = ϵ₃ = @. -0.499ϵ₁
        else
            ϵ₂ = ϵ₃ = Vector{eltype(ϵ₁)}(undef, length(ϵ₁))
        end
        ϵ = collect.(zip(ϵ₁, ϵ₂, ϵ₃))
        σ = collect.(zip(σ₁))
        data = StructArray{BCJMetalDataEntry}((ϵ, σ))
        new{eltype(eltype(ϵ)), eltype(eltype(σ))}(data, name)
    end
    function BCJMetalUniaxialTest(ϵ₁; name, incompressible = true)
        if incompressible
            # ϵ₂ = ϵ₃ = @. sqrt(1 / ϵ₁)
            ϵ₂ = ϵ₃ = @. -0.499ϵ₁
        else
            ϵ₂ = ϵ₃ = Vector{eltype(ϵ₁)}(undef, length(ϵ₁))
        end
        ϵ = collect.(zip(ϵ₁, ϵ₂, ϵ₃))
        σ = collect.(zip(Vector{eltype(ϵ₁)}(undef, length(ϵ₁))))
        data = StructArray{BCJMetalDataEntry}((ϵ, σ))
        new{eltype(eltype(ϵ)), eltype(eltype(σ))}(data, name)
    end
end

include("Bammann1990Modeling.jl")
include("DK.jl")

# end # end of module
