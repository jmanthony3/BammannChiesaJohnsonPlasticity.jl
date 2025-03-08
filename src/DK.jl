using ContinuumMechanicsBase
using ComponentArrays, StructArrays
using Tensors

export DK
export update, predict

struct DK{T<:AbstractFloat} <: AbstractBCJMetalModel
# struct DK{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: AbstractBCJMetalModel
    θ       ::T         # applied temperature
    ϵ̇_eff   ::T         # strain rate (effective)
    ϵₙ      ::T         # final strain
    μ       ::T         # shear modulus at temperature, θ
    N       ::Integer   # number of strain increments
    Δϵ      ::Vector{T} # S         # total strain tensor step
    Δt      ::T         # time step
end

function DK(bcj::BCJMetalStrainControl, μ::AbstractFloat)
    θ       = bcj.θ
    ϵ_dot   = bcj.ϵ_dot
    ϵₙ      = bcj.ϵₙ
    N       = bcj.N
    loadtype= bcj.loadtype
    M       = N + 1
    T       = typeof(float(θ))
    Δϵ      = zeros(T, 6)       # strain increment
    # S       = SymmetricTensor{2, 3, T}
    # Δϵ      = zero(S) # strain increment
    Δt      = (ϵₙ / N) / ϵ_dot

    # state evaluation - loading type
    ϵ̇_eff = if loadtype ∈ (:tension, :compression)    # uniaxial tension/compression
        δϵ  = ϵₙ / N
        Δϵ .= [δϵ, -0.499δϵ, -0.499δϵ, 0., 0., 0.]
        # Δϵ  = S([δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ])
        Δt  = δϵ / ϵ_dot            # timestep
        ϵ_dot
    elseif loadtype == :torsion                                 # torsion
        # convert equivalent strain to true shear strain
        ϵₙ *= 0.5 * √(3.)
        Δϵ .= [0., 0., 0., ϵₙ / N, 0., 0.]
        # Δϵ  = S([0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0])
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[4] / ϵ_dot         # timestep
        # Δt  = Δϵ[1, 2] / ϵ_dot      # timestep
        2ϵ_dot / √3.
    end
    return DK{T}(θ, ϵ̇_eff, ϵₙ, μ, N, Δϵ, Δt)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function update(model::DK, σ__, α__, κ, ϵ__, ϵₚ__, (;
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
    θ       = model.θ
    ϵ̇_eff   = model.ϵ̇_eff
    μ       = model.μ
    ϵₙ      = model.ϵₙ
    N       = model.N
    Δϵ      = model.Δϵ
    Δt      = model.Δt
    M       = N + 1
    T       = typeof(float(θ))
    sqrt23  = √(2.0 / 3.0)
    # if ad_type != AutoForwardDiff()
    #     S       = SymmetricTensor{2, 3} # , T}
    #     σ__     = S(σ__)
    #     α__     = S(α__)
    #     ϵ__     = S(ϵ__)
    #     ϵₚ__    = S(ϵₚ__)
    # end
    # temperature dependent constants
    V   = C₁    * exp( -C₂ / θ )
    Y   = C₃    * exp(  C₄ / θ )
    f   = C₅    * exp( -C₆ / θ )
    β   = Y + (V * asinh( ϵ̇_eff / f ))
    r_d = C₇    * exp( -C₈  / θ )
    h   = C₉    -    (  C₁₀ * θ )
    r_s = C₁₁   * exp( -C₁₂ / θ )
    R_d = C₁₃   * exp( -C₁₄ / θ )
    H   = C₁₅   -    (  C₁₆ * θ )
    R_s = C₁₇   * exp( -C₁₈ / θ )
    Y  *= (C₁₉ < 0.) ? (1.) : (0.5 * ( 1.0 + tanh(max(0., C₁₉ * ( C₂₀ - θ )))))


    # trial guesses
    α_mag       = symmetricmagnitude(α__)
    # α_mag       = norm(α__)
    α_mag      *= sqrt23
    # trial guesses: ISVs (from recovery) and stress
    recovery    = Δt * (r_d * ϵ̇_eff + r_s) * α_mag    # recovery for alpha (kinematic hardening)
    Recovery    = Δt * (R_d * ϵ̇_eff + R_s) * κ  # recovery for kappa (isotropic hardening)
    αₜᵣ__       = α__   .* (1 - recovery)
    # αₜᵣ__       = α__   * (1 - recovery)
    κₜᵣ         = κ     * (1 - Recovery)
    σₜᵣ__       = σ__ + (2μ * Δϵ)                         # trial stress
    ξ__         = σₜᵣ__ - (2.0 / 3.0) * αₜᵣ__                                             # trial overstress original
    # ξ__         = σₜᵣ__ - αₜᵣ__                                             # trial overstress original
    # ξ__          .= σₜᵣ__ - sqrt23 .* αₜᵣ__                                 # trial overstress FIT
    ξ_mag       = symmetricmagnitude(ξ__)
    # ξ_mag       = norm(ξ__)


    # yield criterion
    flow_rule = ξ_mag - sqrt23 * (κₜᵣ + β)                                             # same as vumat20
    if flow_rule <= 0.      # elastic
        # trial guesses are correct
        σ__             = @. σₜᵣ__
        α__             = @. αₜᵣ__
        κ               = κₜᵣ
        # state.ξ__               = ξ__
        ϵ__            += @. Δϵ
        # state.ϵ_dot_plastic__    .= 0.
    else                    # plastic
        # Radial Return
        Δγ              = flow_rule / (2μ + 2(h + H) / 3)     # original
        n̂               = ξ__ ./ ξ_mag
        # n̂               = ξ__ / ξ_mag
        σ__prev         = σ__
        σ__             = @. σₜᵣ__ - (2μ * Δγ) .* n̂
        α__             = @. αₜᵣ__ + ( h * Δγ) .* n̂
        # σ__           = @. σₜᵣ__ - (2μ * Δγ) * n̂
        # α__           = @. αₜᵣ__ + ( h * Δγ) * n̂
        # state.ξ__       = state.σ__ - state.α__
        κ               = κₜᵣ   + (H * sqrt23 * Δγ)  # original
        ϵₚ__           += @. (Δϵ - ((σ__ - σ__prev) ./ 2μ))
        # ϵₚ__           += @. (Δϵ - ((σ__ - σ__prev) / 2μ))
        ϵ__            += @. Δϵ
    end
    # ϵ_dot_plastic__ .= (f * sinh(V \ (ξ_mag - κ - Y)) / ξ_mag) * ξ__
    return σ__, α__, κ, ϵ__, ϵₚ__
    # return triu_vec(σ__), triu_vec(α__), κ, triu_vec(ϵ__), triu_vec(ϵₚ__)
    # return nothing
end

function ContinuumMechanicsBase.predict(
            ψ   ::DK{T}, # , S},
            test::AbstractBCJMetalTest{T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat} # , S<:SymmetricTensor{2, 3}}
    M = ψ.N + 1
    σ__     = zeros(T, 6)   # deviatoric stress
    ϵₚ__    = zeros(T, 6)   # plastic strain
    ϵ__     = zeros(T, 6)   # total strain
    α__     = fill(1e-7, 6) # alpha: kinematic hardening
    κ       = 0.0           # kappa: isotropic hardening
    ξ__     = zeros(T, 6)   # overstress (S - 2/3*alpha)
    ϵ⃗ = []
    σ⃗ = []
    push!(ϵ⃗, ϵ__)
    push!(σ⃗, σ__)
    for i ∈ range(2, M)
        σ__, α__, κ, ϵ__, ϵₚ__ = update(ψ, σ__, α__, κ, ϵ__, ϵₚ__, p)
        # update!(ψ, σ__, α__, κ, ϵ__, ϵₚ__, p)
        push!(ϵ⃗, ϵ__)
        push!(σ⃗, σ__)
    end
    return (data=(ϵ=hcat(ϵ⃗...), σ=hcat(σ⃗...)),)
    # σ__     = zeros(T, 6)   # deviatoric stress
    # ϵₚ__    = zeros(T, 6)   # plastic strain
    # ϵ__     = zeros(T, 6)   # total strain
    # α__     = fill(1e-7, 6) # alpha: kinematic hardening
    # κ       = 0.0           # kappa: isotropic hardening
    # ξ__     = zeros(T, 6)   # overstress (S - 2/3*alpha)
    # ϵ⃗ = zeros(T, (6, M))
    # σ⃗ = zeros(T, (6, M))
    # for i ∈ range(2, M)
    #     σ__, α__, κ, ϵ__, ϵₚ__ = update(ψ, σ__, α__, κ, ϵ__, ϵₚ__, p)
    #     ϵ⃗[:, i], σ⃗[:, i] = ϵ__, σ__
    # end
    # s = SymmetricTensor{2, 3, T}
    # if ad_type != AutoFiniteDiff()
    #     σ__     = zero(s)       # deviatoric stress
    #     ϵₚ__    = zero(s)       # plastic strain
    #     ϵ__     = zero(s)       # total strain
    #     α__     = fill(1e-7, s) # alpha: kinematic hardening
    #     κ       = 0.            # kappa: isotropic hardening
    #     ξ__     = zero(s)       # overstress (S - 2/3*alpha)
    #     ϵ⃗ = zeros(s, M, 1)
    #     σ⃗ = zeros(s, M, 1)
    #     ϵ⃗[1], σ⃗[1] = ϵ__, σ__
    #     for i ∈ range(2, M)
    #         σ__, α__, κ, ϵ__, ϵₚ__ = update(ψ, σ__, α__, κ, ϵ__, ϵₚ__, p)
    #         ϵ⃗[i], σ⃗[i] = s(ϵ__), s(σ__)
    #     end
    # else
    # end
    # return (data=(ϵ=ϵ⃗, σ=σ⃗),)
end