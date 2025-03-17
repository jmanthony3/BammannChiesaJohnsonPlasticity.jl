module Metals



export AbstractBCJMetalModel
export BCJMetalStrainControl
export BCJMetalDataEntry
export BCJMetalUniaxialTest

import ..BammannChiesaJohnsonPlasticity: AbstractBCJModel, AbstractBCJTest

using ContinuumMechanicsBase
using ComponentArrays, StructArrays
using DocStringExtensions
# using LinearAlgebra


"Sub-type for BCJ-models specific to metals."
abstract type AbstractBCJMetalModel         <: AbstractBCJModel end
# abstract type ISVMetal{T<:BCJMetal} end
# abstract type ISVMetalKinematicHardening{T} <: ISVMetal{T} end # α__
# abstract type ISVMetalIsotropicHardening{T} <: ISVMetal{T} end # κ
# abstract type ISVMetalDamage{T}             <: ISVMetal{T} end # ϕ
"Sub-type for BCJ-models specific to metals."
abstract type AbstractBCJMetalTest{T}       <: AbstractBCJTest end

"Stucture for strain-controlled loadings of metals for temperature, `θ`; strain rate, `ϵ̇`; final strain, `ϵₙ`; number of loading increments, `N`; loading direction, `loaddir` ∈ {`:tension`, `:compression`, `:torsion`}."
struct BCJMetalStrainControl{T1<:Integer, T2<:AbstractFloat} <: AbstractBCJTest
    θ       ::T2        # applied temperature
    ϵ̇       ::T2        # applied strain rate
    ϵₙ      ::T2        # final strain
    N       ::T1        # number of strain increments
    loaddir ::Symbol    # load type (:tension, :compression, :torsion)
end

"Store vectors of strain, `ϵ`, and stress, `σ` data."
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

end # end of module
