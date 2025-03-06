# using PlasticityBase
using ContinuumMechanicsBase
using Tensors

export DK, update

mutable struct DK{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: ContinuumMechanicsBase.AbstractMaterialModel
    θ               ::T         # applied temperature
    μ               ::T         # shear modulus at temperature, θ
    σ__             ::Vector{T} # deviatoric stress tensor
    ϵₚ__            ::Vector{T} # plastic strain tensor
    ϵ_dot_plastic__ ::Vector{T} # plastic strain rate
    ϵ__             ::Vector{T} # total strain tensor
    ϵ_dot_effective ::T         # strain rate (effective)
    Δϵ              ::Vector{T} # total strain tensor step
    Δt              ::T         # time step
    α__             ::Vector{T} # kinematic hardening tensor
    κ               ::T         # isotropic hardening scalar
    ξ__             ::Vector{T} # overstress tensor (S - 2/3*alpha)
end

function DK(bcj::BCJMetalStrainControl, μ::AbstractFloat)
    θ       = bcj.θ
    ϵ_dot   = bcj.ϵ_dot
    ϵₙ      = bcj.ϵₙ
    N       = bcj.N
    loadtype= bcj.loadtype
    M       = N + 1
    T       = typeof(float(θ))
    S       = SymmetricTensor{2, 3, T}
    # array declarations
    ## OSVs
    σ__             = zero(S)       # deviatoric stress
    ϵₚ__            = zero(S)       # plastic strain
    ϵ_dot_plastic__ = zero(S)       # plastic strain rate
    ϵ__             = zero(S)       # total strain
    ## ISVs
    α__             = fill(1e-7, S) # alpha: kinematic hardening
    κ               = 0.            # kappa: isotropic hardening
    ## holding values
    Δϵ              = zero(S)       # strain increment
    ξ__             = zero(S)       # overstress (S - 2/3*alpha)


    # state evaluation - loading type
    ϵ_dot_effective = if loadtype ∈ (:tension, :compression)    # uniaxial tension/compression
        δϵ  = ϵₙ / N
        Δϵ  = S([δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ])
        Δt  = δϵ / ϵ_dot            # timestep
        ϵ_dot
    elseif loadtype == :torsion                                 # torsion
        # convert equivalent strain to true shear strain
        ϵₙ *= 0.5 * √(3.)
        Δϵ  = S([0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0])
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[1, 2] / ϵ_dot      # timestep
        2ϵ_dot / √3.
    end
    return DK{T}(N, θ, μ,
        σ__, ϵₚ__, ϵ_dot_plastic__, ϵ__, ϵ_dot_effective, Δϵ, Δt,
        V, Y, f, h, r_d, r_s, H, R_d, R_s, α__, κ, β, ξ__, σₜᵣ__, αₜᵣ__, κₜᵣ)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function ContinuumMechanicsBase.predict(state::DK, bcj::BCJMetalStrainControl, (
            C₁,     C₂,     # V
            C₃,     C₄,     # Y
            C₅,     C₆,     # f
            C₇,     C₈,     # r_d
            C₉,     C₁₀,    # h
            C₁₁,    C₁₂,    # r_s
            C₁₃,    C₁₄,    # R_d
            C₁₅,    C₁₆,    # H
            C₁₇,    C₁₈,    # R_s
            C₁₉,    C₂₀     # Y_adj
        ))
    θ, μ, Δϵ, Δt    = state.θ, state.μ, state.Δϵ, state.Δt
    ϵ_dot_effective = state.ϵ_dot_effective
    # V, Y, f, β      = bcj.V, bcj.Y, bcj.f, bcj.β
    # h, r_d, r_s     = bcj.h, bcj.r_d, bcj.r_s # alpha
    # H, R_d, R_s     = bcj.H, bcj.R_d, bcj.R_s # kappa
    # temperature dependent constants
    V   = C1    * exp( -C2 / θ )
    Y   = C3    * exp(  C4 / θ )
    f   = C5    * exp( -C6 / θ )

    β   = Y + (V * asinh( ϵ_dot_effective / f ))

    r_d = C7    * exp( -C8  / θ )
    h   = C9    -    (  C10 * θ )
    r_s = C11   * exp( -C12 / θ )

    R_d = C13   * exp( -C14 / θ )
    H   = C15   -    (  C16 * θ )
    R_s = C17   * exp( -C18 / θ )

    Y  *= (C19 < 0.) ? (1.) : (0.5 * ( 1.0 + tanh(max(0., C19 * ( C20 - θ )))))

    sqrt23          = √(2 / 3)


    # timestep calculations
    for i ∈ range(2, bcj.N + 1)
        α_mag       = norm(state.α__)
        # α_mag = sqrt( α_mag * 3./2.)       # match cho
        α_mag      *= sqrt23       # match vumat20
        # trial guesses: ISVs (from recovery) and stress
        recovery    = Δt * (r_d * ϵ_dot_effective + r_s) * α_mag    # recovery for alpha (kinematic hardening)
        Recovery    = Δt * (R_d * ϵ_dot_effective + R_s) * state.κ    # recovery for kappa (isotropic hardening)
        αₜᵣ__       = state.α__   * (1 - recovery)
        κₜᵣ         = state.κ     * (1 - Recovery)

        ## trial stress guess
        σₜᵣ__       = state.σ__ + 2μ * Δϵ         # trial stress
        ξ__         = σₜᵣ__ - (2. / 3.) * αₜᵣ__ # trial overstress original
        # ξ__          .= σₜᵣ__ - sqrt23 .* αₜᵣ__ # trial overstress FIT
        ξ_mag       = norm(state.ξ__)



        # ----------------------------------- #
        ###   ---   YIELD CRITERION   ---   ###
        # ----------------------------------- #
        flow_rule = ξ_mag - sqrt23 * (κₜᵣ + β)         # same as vumat20
        # Crit = Xi_mag - (Katr + β) #changed to FIT
        if flow_rule <= 0.      # elastic
            # trial guesses are correct
            state.σ__                 = σₜᵣ__
            state.α__                 = αₜᵣ__
            state.κ                   = κₜᵣ
            state.ξ__               = ξ__
            state.ϵ__                += Δϵ
            # state.ϵ_dot_plastic__    .= 0.
        else                    # plastic
            # Radial Return
            Δγ                      = flow_rule / (2μ + 2(h + H) / 3)     # original
            n̂                       = state.ξ__ / ξ_mag
            σ__prev                 = state.σ__
            state.σ__               = σₜᵣ__ - (2μ * Δγ) * n̂
            state.α__               = αₜᵣ__ + ( h * Δγ) * n̂
            state.ξ__               = state.σ__ - state.α__
            state.κ                 = κₜᵣ   + (H * sqrt23 * Δγ)  # original
            state.ϵₚ__             += (Δϵ - ((state.σ__ - σ__prev) / 2μ))
            state.ϵ__              += Δϵ
        end
        state.ϵ_dot_plastic__    .= (f * sinh(V \ (ξ_mag - state.κ - Y)) / ξ_mag) .* state.ξ__
        # record!(history, i, bcj)
    end
    return nothing
end