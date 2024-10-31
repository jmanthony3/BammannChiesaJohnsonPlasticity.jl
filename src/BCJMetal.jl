abstract type BCJMetal              <: BCJ end
abstract type ISVMetal{T<:BCJMetal} end
abstract type KinematicHardening{T} <: ISVMetal{T} end # α__
abstract type IsotropicHardening{T} <: ISVMetal{T} end # κ
abstract type Damage{T}             <: ISVMetal{T} end # ϕ

symmetricmagnitude(tensor::Vector{<:Real}) = √( sum(tensor[1:3] .^ 2.) + 2sum(tensor[4:6] .^ 2.) )

function symmetricvonMises(tensor::Union{Vector{<:Real}, SubArray{<:Real}})::AbstractFloat
    σvM = sum(map(x->x^2., [tensor[1] - tensor[2], tensor[2] - tensor[3], tensor[3] - tensor[1]])) + (
        6sum(map(x->x^2., [tensor[4], tensor[5], tensor[6]])))
    return √(σvM / 2.)
end
function symmetricvonMises(tensor::Matrix{<:Real})::Vector{AbstractFloat}
    return map(symmetricvonMises, eachcol(tensor))
end

struct BCJMetalStrainControl{T1<:Integer, T2<:AbstractFloat} <: BCJMetal
    θ       ::T2                # applied temperature
    ϵ_dot   ::T2                # applied strain rate
    ϵₙ      ::T2                # final strain
    N       ::T1                # number of strain increments
    istate  ::T1                # load type (1: uniaxial tension; 2: torsion)
    params  ::Dict{String, T2}  # material constants
end

mutable struct BCJMetalCurrentConfiguration{Version<:BCJMetal, T<:AbstractFloat} <: BCJMetal
    N               ::Integer   # number of strain increments
    θ               ::T         # applied temperature
    μ               ::T         # shear modulus at temperature, θ
    σ__             ::Vector{T} # deviatoric stress tensor
    ϵₚ__            ::Vector{T} # plastic strain tensor
    ϵ_dot_plastic__ ::Vector{T} # plastic strain rate
    ϵ__             ::Vector{T} # total strain tensor
    ϵ_dot_effective ::T         # strain rate (effective)
    Δϵ              ::Vector{T} # total strain tensor step
    Δt              ::T         # time step
    V               ::T         # strain rate sensitivity of yield stress at temperature, θ
    Y               ::T         # rate independent yield stress at temperature, θ
    f               ::T         # strain rate at which yield becomes strain rate dependent at temperature, θ
    h               ::T         # kinematic hardening modulus at temperature, θ
    r_d             ::T         # dynamic recovery of kinematic hardening at temperature, θ
    r_s             ::T         # diffusion controlled static/thermal recovery of kinematic hardening at temperature, θ
    H               ::T         # isotropic hardening modulus at temperature, θ
    R_d             ::T         # dynamic recovery of isotropic hardening at temperature, θ
    R_s             ::T         # diffusion controlled static/thermal recovery of isotropic hardening at temperature, θ
    α__             ::Vector{T} # kinematic hardening tensor
    κ               ::T         # isotropic hardening scalar
    β               ::T         # yield function
    ξ__             ::Vector{T} # overstress tensor (S - 2/3*alpha)
    σₜᵣ__           ::Vector{T} # deviatoric stress tensor (trial)
    αₜᵣ__           ::Vector{T} # kinematic hardening tensor (trial)
    κₜᵣ             ::T         # isotropic hardening (trial)
end

mutable struct BCJMetalConfigurationHistory{T<:AbstractFloat} <: BCJMetal
    σ__             ::Matrix{T} # deviatoric stress tensor
    ϵₚ__            ::Matrix{T} # plastic strain tensor
    ϵ_dot_plastic__ ::Matrix{T} # plastic strain rate
    ϵ__             ::Matrix{T} # total strain tensor
    α__             ::Matrix{T} # kinematic hardening tensor
    κ               ::Vector{T} # isotropic hardening scalar
    ξ__             ::Matrix{T} # overstress tensor (S - 2/3*alpha)
end

function Base.:+(x::T, y::T) where {T<:BCJMetalConfigurationHistory}
    return BCJMetalConfigurationHistory{eltype(x.σ__)}(
        hcat(x.σ__,                y.σ__),
        hcat(x.ϵₚ__,               y.ϵₚ__),
        hcat(x.ϵ_dot_plastic__,    y.ϵ_dot_plastic__),
        hcat(x.ϵ__,                y.ϵ__ .+ x.ϵ__[:, end]),
        hcat(x.α__,                y.α__),
        vcat(x.κ,                  y.κ),
        hcat(x.ξ__,                y.ξ__)
    )
end

function Base.copyto!(reference::BCJMetalCurrentConfiguration, history::BCJMetalConfigurationHistory)
    # for attr ∈ (:θ, :V, :Y, :f, :h, :r_d, :r_s, :H, :R_d, :R_s, :α__, :κ, :β, :ξ__)
    #     setfield!(reference, attr, getfield(current, attr))
    # end
    reference.σ__               = history.σ__[:, end]
    # reference.ϵₚ__              = history.ϵₚ__[:, end]
    # reference.ϵ_dot_plastic__   = history.ϵ_dot_plastic__[:, end]
    # reference.ϵ__               = history.ϵ__[:, end]
    reference.α__               = history.α__[:, end]
    reference.κ                 = history.κ[end]
    reference.ξ__               = history.ξ__[:, end]
    return nothing
end

function record!(history::BCJMetalConfigurationHistory, i::Integer, current::BCJMetalCurrentConfiguration)
    history.σ__[:, i]              .= current.σ__
    history.ϵₚ__[:, i]             .= current.ϵₚ__
    history.ϵ_dot_plastic__[:, i]  .= current.ϵ_dot_plastic__
    history.ϵ__[:, i]              .= current.ϵ__
    history.α__[:, i]              .= current.α__
    history.κ[i]                    = current.κ
    history.ξ__[:, i]              .= current.ξ__
    return nothing
end