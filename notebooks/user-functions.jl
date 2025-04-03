# module B93F

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays, StructArrays
# using Tensors # uncomment when we can work with Tensors.jl
using DocStringExtensions

# """
# Structure for viscoplasticity model with loading conditions and material properties.
# Here, uses the effective strain rate based on applied strain rate and loading direction.
# """
# struct Bammann1993Failure{T<:AbstractFloat} <: BammannChiesaJohnsonPlasticity.AbstractBCJMetalModel
# # struct Bammann1993Failure{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: AbstractBCJMetalModel
#     θ       ::T         # applied temperature
#     ϵ̇_eff   ::T         # strain rate (effective)
#     ϵₙ      ::T         # final strain
#     N       ::Integer   # number of strain increments
#     Δϵ      ::Vector{T} # S         # total strain tensor step
#     Δt      ::T         # time step
#     μ       ::T         # shear modulus at temperature, θ
# end

"""
    $(SIGNATURES)

Outer constructor for loading conditions and material properties which assumes a Poisson's ratio of 0.5.
Here, `μ` is the shear modulus.
"""
function Bammann1993Failure(Ω::BammannChiesaJohnsonPlasticity.BCJMetalStrainControl, μ::AbstractFloat)
    θ       = Ω.θ
    ϵ̇       = Ω.ϵ̇
    ϵₙ      = Ω.ϵₙ
    N       = Ω.N
    loaddir = Ω.loaddir
    M       = N + 1
    T       = typeof(float(θ))
    Δϵ      = zeros(T, 6)       # strain increment
    # S       = SymmetricTensor{2, 3, T}
    # Δϵ      = zero(S) # strain increment
    Δt      = (ϵₙ / N) / ϵ̇

    # state evaluation - loading type
    ϵ̇_eff = if loaddir ∈ (:tension, :compression)    # uniaxial tension/compression
        δϵ  = ϵₙ / N
        Δϵ .= [δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ]
        # Δϵ  = S([δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ])
        Δt  = δϵ / ϵ̇            # timestep
        ϵ̇
    elseif loaddir == :torsion                                 # torsion
        # convert equivalent strain to true shear strain
        ϵₙ *= 0.5 * √(3.)
        Δϵ .= [0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0]
        # Δϵ  = S([0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0])
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[2] / ϵ̇         # timestep
        # Δt  = Δϵ[1, 2] / ϵ_dot      # timestep
        2ϵ̇ / √3.
    end
    return Bammann1993Failure{T}(θ, ϵ̇_eff, ϵₙ, N, Δϵ, Δt, μ)
end

"""
Using the equations and constants from [Bammann et. al. (1993)](@cite bammannFailureDuctileMaterials1993), this kernel function maps the current material state and ISVs onto the next configuration.
Though not explicitly listed in paper, temperature equations `h = C₁₅ * exp(-C₁₆ / θ)` and `H = C₁₇ * exp(-C₁₈ / θ)` are included (and their constants renumbered) from (c. f. [Horstemeyer (1994)](@cite horstemeyerPredictingFormingLimit1994)).
Important: `ϕ` is included in the list of arguments, but is presently, internally set to zero.
This is a limitation of the point simulator causing infinite stress triaxiality, χ.
"""
function update(ψ::Bammann1993Failure, σ̲̲, α̲̲, κ, ϕ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, (;
            C₁,     C₂,     # V
            C₃,     C₄,     # Y
            C₅,     C₆,     # f
            C₇,     C₈,     # r_d
            C₉,     C₁₀,    # r_s
            C₁₁,    C₁₂,    # R_d
            C₁₃,    C₁₄,    # R_s
            C₁₅,    C₁₆,    # h
            C₁₇,    C₁₈,    # H
            m̄               # ϕ
        ))
    θ       = ψ.θ
    ϵ̇_eff   = ψ.ϵ̇_eff
    μ       = ψ.μ
    ϵₙ      = ψ.ϵₙ
    N       = ψ.N
    Δϵ      = ψ.Δϵ
    Δt      = ψ.Δt
    M       = N + 1
    T       = typeof(float(θ))
    sqrt23  = √(2.0 / 3.0)
    ϕ       = 0.0 # enforce zero damage for point simulator
    ϕ̇       = 0.0
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
    # β   = ( Y * (1 - ϕ) ) + (  V * (1 - ϕ) * asinh( ϵ̇_eff / f )  )
    # (^) will need this eventually | (v) this is for point simulator
    β   = Y + (V * asinh( ϵ̇_eff / f ))
    # @show V, Y, f, ϕ, β
    r_d = C₇    * exp( -C₈  / θ )
    r_s = C₉    * exp( -C₁₀ / θ )
    R_d = C₁₁   * exp( -C₁₂ / θ )
    R_s = C₁₃   * exp( -C₁₄ / θ )
    h   = C₁₅   * exp( -C₁₆ / θ )
    H   = C₁₇   * exp( -C₁₈ / θ )
    # @show r_d, r_s, R_d, R_s, h, H
    p(σ)= sum(σ[[1, 4, 6]]) / 3.0
    # # if p(σ̲̲) == 0.0
    # #     χ, ϕ̇ = 0.0, 0.0
    # # else
    # # end
    # # num = ( 2(2m̄ - 1) * p(σ̲̲) )
    # # den = ( (2m̄ + 1) * vonMises(σ̲̲) )
    # χ   = sinh(#=[=#  ( 2(2m̄ - 1) * p(σ̲̲) ) / ( (2m̄ + 1) * vonMises(σ̲̲) )  #=]=#)
    # ϕ̇   = χ * (#=[=#(   1  /  ( (1 - ϕ)^m̄ )   )   -   (1 - ϕ)#=]=#) * ϵ̇_eff
    # # @show num, den, χ, ϕ̇


    # trial guesses
    α_mag       = norm_symvec(α̲̲)
    α_mag      *= sqrt23
    # α_mag       = norm(α̲̲)
    # trial guesses: ISVs (from recovery) and stress
    recovery    = Δt * (r_d * ( sqrt23 * ϵ̇_eff ) + r_s) * α_mag     # recovery for alpha (kinematic hardening)
    Recovery    = Δt * (R_d * ( sqrt23 * ϵ̇_eff ) + R_s) * κ         # recovery for kappa (isotropic hardening)
    # σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ .* (1 - ϕ) .* Δϵ) .- ( (ϕ̇ * Δt) / (1 - ϕ) )  # deviatoric stress (trial)
    # (^) will need this eventually | (v) this is for point simulator
    σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ .* Δϵ)                                    # deviatoric stress (trial)
    σ̲̲′⁽ᵗʳ⁾      = σ̲̲⁽ᵗʳ⁾ - (p(σ̲̲⁽ᵗʳ⁾) .* [1, 0, 0, 1, 0, 1])
    α̲̲⁽ᵗʳ⁾       = α̲̲    .* (1 - recovery)
    # α̲̲⁽ᵗʳ⁾       = α̲̲     * (1 - recovery)
    κ⁽ᵗʳ⁾       = κ     * (1 - Recovery)
    ξ̲̲⁽ᵗʳ⁾       = σ̲̲′⁽ᵗʳ⁾ - (2.0 / 3.0 * α̲̲⁽ᵗʳ⁾)                      # over-stress (trial)
    ξ_mag       = norm_symvec(ξ̲̲⁽ᵗʳ⁾)
    # ξ_mag       = norm(ξ__)


    # yield criterion
    F = ξ_mag - sqrt23 * (κ⁽ᵗʳ⁾ + β) # * (1 - ϕ)
    # @show F
    if F <= 0.0     # elastic
        # trial guesses are correct
        σ̲̲       = @. σ̲̲⁽ᵗʳ⁾
        α̲̲       = @. α̲̲⁽ᵗʳ⁾
        κ       =    κ⁽ᵗʳ⁾
        ϕ       =    ϕ
        ϵ̲̲      += @. Δϵ
        # state.ϵ_dot_plastic__    .= 0.
    else            # plastic
        # Radial Return
        Δγ      = F / (2μ + 3\2(h + H))
        n̂       = ξ̲̲⁽ᵗʳ⁾ ./ ξ_mag
        # n̂       = ξ__ / ξ_mag
        σ̲̲_prev  = σ̲̲
        σ̲̲       = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) .* n̂
        α̲̲       = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) .* n̂
        # σ̲̲       = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) * n̂
        # α̲̲       = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) * n̂
        # κ       =    κ⁽ᵗʳ⁾ + ( H * Δγ * sqrt23)
        # χ       = sinh(#=[=#  ( 2(2m̄ - 1) * p(σ̲̲) ) / ( (2m̄ + 1) * vonMises(σ̲̲) )  #=]=#)
        # ϕ       =    1 - (#={=#
        #         1 + (#=[=#   (1 - ϕ) ^ (1 + m̄) - 1   #=]=#) * exp(#=[=#
        #             (#=d̄: =# sqrt23 * ϵ̇_eff) * χ * (1 + m̄) * Δt   #=]=#)
        #     #=}=#) ^ ( 1 / (1 + m̄) )
        # # (^) exact solution | (v) Forward-Euler method
        # # ϕ̇       = χ * (#=[=#(   1  /  ( (1 - ϕ)^m̄ )   )   -   (1 - ϕ)#=]=#) * ϵ̇_eff
        # # ϕ      +=    ϕ̇ * Δt
        # @show p(σ̲̲), vonMises(σ̲̲), χ, ϕ
        ϕ       = ϕ
        ϵ̲̲⁽ᵖ⁾   += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) ./ 2μ))
        # ϵ̲̲⁽ᵖ⁾   += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) / 2μ))
        ϵ̲̲      += @. Δϵ
    end
    # ϵ_dot_plastic__  = @. (f * sinh(V \ (ξ_mag - κ - Y)) / ξ_mag) * ξ__
    if ϕ > 0.99
        error("Element death.")
    end
    return σ̲̲, α̲̲, κ, ϕ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾
    # return triu_vec(σ__), triu_vec(α__), κ, triu_vec(ϵ__), triu_vec(ϵₚ__)
    # return nothing
end

function ContinuumMechanicsBase.predict(
            ψ   ::Bammann1993Failure{T}, # , S},
            test::BammannChiesaJohnsonPlasticity.AbstractBCJMetalTest{T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat} # , S<:SymmetricTensor{2, 3, T}}
    M = ψ.N + 1
    σ̲̲       = zeros(T, 6)   # deviatoric stress
    ϵ̲̲⁽ᵖ⁾    = zeros(T, 6)   # plastic strain
    ϵ̲̲       = zeros(T, 6)   # total strain
    α̲̲       = fill(1e-7, 6) # kinematic hardening
    κ       = 0.0           # isotropic hardening
    ϕ       = 0.0           # damage
    ϵ⃗ = []
    σ⃗ = []
    push!(ϵ⃗, ϵ̲̲)
    push!(σ⃗, σ̲̲)
    for i ∈ range(2, M)
        σ̲̲, α̲̲, κ, ϕ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾ = update(ψ, σ̲̲, α̲̲, κ, ϕ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, p)
        # update!(ψ, σ__, α__, κ, ϵ__, ϵₚ__, p)
        push!(ϵ⃗, ϵ̲̲)
        push!(σ⃗, σ̲̲)
    end
    return (data=(ϵ=hcat(ϵ⃗...), σ=hcat(σ⃗...)),)
    # σ__     = zeros(T, 6)   # deviatoric stress
    # ϵₚ__    = zeros(T, 6)   # plastic strain
    # ϵ__     = zeros(T, 6)   # total strain
    # α__     = fill(1e-7, 6) # alpha: kinematic hardening
    # κ       = 0.0           # kappa: isotropic hardening
    # ϕ       = 0.0           # phi: damage
    # ϵ⃗ = zeros(T, (6, M))
    # σ⃗ = zeros(T, (6, M))
    # for i ∈ range(2, M)
    #     σ__, α__, κ, ϕ, ϵ__, ϵₚ__ = update(ψ, σ__, α__, κ, ϕ, ϵ__, ϵₚ__, p)
    #     ϵ⃗[:, i], σ⃗[:, i] = ϵ__, σ__
    # end
    # s = SymmetricTensor{2, 3, T}
    # if ad_type != AutoFiniteDiff()
    #     σ__     = zero(s)       # deviatoric stress
    #     ϵₚ__    = zero(s)       # plastic strain
    #     ϵ__     = zero(s)       # total strain
    #     α__     = fill(1e-7, s) # alpha: kinematic hardening
    #     κ       = 0.            # kappa: isotropic hardening
    #     ϕ       = 0.            # phi: damage
    #     ϵ⃗ = zeros(s, M, 1)
    #     σ⃗ = zeros(s, M, 1)
    #     ϵ⃗[1], σ⃗[1] = ϵ__, σ__
    #     for i ∈ range(2, M)
    #         σ__, α__, κ, ϕ, ϵ__, ϵₚ__ = update(ψ, σ__, α__, κ, ϕ, ϵ__, ϵₚ__, p)
    #         ϵ⃗[i], σ⃗[i] = s(ϵ__), s(σ__)
    #     end
    # else
    # end
    # return (data=(ϵ=ϵ⃗, σ=σ⃗),)
end

"""
Constants for temperature equations from [Bammann et. al. (1993)](@cite bammannFailureDuctileMaterials1993).
Note: though not explicitly listed in paper, temperature equations `h = C₁₅ * exp(-C₁₆ / θ)` and `H = C₁₇ * exp(-C₁₈ / θ)` are included (and their constants renumbered) from (c. f. [Horstemeyer (1994)](@cite horstemeyerPredictingFormingLimit1994)).
"""
ContinuumMechanicsBase.parameters(::Bammann1993Failure) = (
    :C₁,    :C₂,    # V
    :C₃,    :C₄,    # Y
    :C₅,    :C₆,    # f
    :C₇,    :C₈,    # r_d
    :C₉,    :C₁₀,   # r_s
    :C₁₁,   :C₁₂,   # R_d
    :C₁₃,   :C₁₄,   # R_s
    :C₁₅,   :C₁₆,   # h
    :C₁₇,   :C₁₈,   # H
    :m̄              # ϕ
)

# end # end of module
