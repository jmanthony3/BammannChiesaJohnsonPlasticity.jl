export Bammann1990Modeling
export update, predict

using ContinuumMechanicsBase
using ComponentArrays, StructArrays
# using Tensors # uncomment when we can work with Tensors.jl

"""
Structure for viscoplasticity model with loading conditions and material properties.
Here, uses the effective strain rate based on applied strain rate and loading direction.
"""
struct Bammann1990Modeling{T<:AbstractFloat} <: AbstractBCJMetalModel
# struct Bammann1990Modeling{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: AbstractBCJMetalModel
    θ       ::T         # applied temperature
    ϵ̇_eff   ::T         # strain rate (effective)
    ϵₙ      ::T         # final strain
    N       ::Integer   # number of strain increments
    Δϵ      ::Vector{T} # S         # total strain tensor step
    Δt      ::T         # time step
    μ       ::T         # shear modulus at temperature, θ
end

"""
    $(SIGNATURES)

Outer constructor for loading conditions and material properties which assumes a Poisson's ratio of 0.5.
Here, `μ` is the shear modulus.
"""
function Bammann1990Modeling(Ω::BCJMetalStrainControl, μ::AbstractFloat)
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
    return Bammann1990Modeling{T}(θ, ϵ̇_eff, ϵₙ, N, Δϵ, Δt, μ)
end

"""
Using the equations and constants from [Bammann (1990)](@cite bammannModelingTemperatureStrain1990), this kernel function maps the current material state and ISVs onto the next configuration.
Note: though not explicitly listed in paper, temperature equations `h = C₁₅ * exp(C₁₆ / θ)` and `H = C₁₇ * exp(C₁₈ / θ)` are included (c. f. [DYNA3D User Manual (1993)](@cite whirley1993dyna3d)).
"""
function update(ψ::Bammann1990Modeling, σ̲̲, α̲̲, κ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, (;
            C₁,     C₂,     # V
            C₃,     C₄,     # Y
            C₅,     C₆,     # f
            C₇,     C₈,     # r_d
            C₉,     C₁₀,    # r_s
            C₁₁,    C₁₂,    # R_d
            C₁₃,    C₁₄,    # R_s
            C₁₅,    C₁₆,    # h
            C₁₇,    C₁₈     # H
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
    r_s = C₉    * exp( -C₁₀ / θ )
    R_d = C₁₁   * exp( -C₁₂ / θ )
    R_s = C₁₃   * exp( -C₁₄ / θ )
    h   = C₁₅   * exp(  C₁₆ / θ )
    H   = C₁₇   * exp(  C₁₈ / θ )


    # trial guesses
    α_mag       = norm_symvec(α̲̲)
    # α_mag       = norm(α̲̲)
    # trial guesses: ISVs (from recovery) and stress
    recovery    = Δt * (r_d * ϵ̇_eff + r_s) * α_mag      # recovery for alpha (kinematic hardening)
    Recovery    = Δt * (R_d * ϵ̇_eff + R_s) * κ          # recovery for kappa (isotropic hardening)
    α̲̲⁽ᵗʳ⁾       = α̲̲    .* (1 - recovery)
    # α̲̲⁽ᵗʳ⁾       = α̲̲     * (1 - recovery)
    κ⁽ᵗʳ⁾       = κ     * (1 - Recovery)
    σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ * Δϵ)                         # deviatoric stress (trial)
    ξ̲̲⁽ᵗʳ⁾       = σ̲̲⁽ᵗʳ⁾ - α̲̲⁽ᵗʳ⁾                         # over-stress (trial)
    ξ_mag       = norm_symvec(ξ̲̲⁽ᵗʳ⁾)
    # ξ_mag       = norm(ξ__)


    # yield criterion
    F = ξ_mag - κ⁽ᵗʳ⁾ - β
    if F <= 0.0     # elastic
        # trial guesses are correct
        σ̲̲       = @. σ̲̲⁽ᵗʳ⁾
        α̲̲       = @. α̲̲⁽ᵗʳ⁾
        κ       = κ⁽ᵗʳ⁾
        ϵ̲̲      += @. Δϵ
        # state.ϵ_dot_plastic__    .= 0.
    else            # plastic
        # Radial Return
        Δγ      = F / (2μ + 2(h + H) / 3)
        n̂       = ξ̲̲⁽ᵗʳ⁾ ./ ξ_mag
        # n̂       = ξ__ / ξ_mag
        σ̲̲_prev  = σ̲̲
        σ̲̲       = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) .* n̂
        α̲̲       = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) .* n̂
        # σ̲̲       = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) * n̂
        # α̲̲       = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) * n̂
        κ       =    κ⁽ᵗʳ⁾ + ( H * Δγ)
        ϵ̲̲⁽ᵖ⁾   += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) ./ 2μ))
        # ϵₚ__   += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) / 2μ))
        ϵ̲̲      += @. Δϵ
    end
    # ϵ_dot_plastic__  = @. (f * sinh(V \ (ξ_mag - κ - Y)) / ξ_mag) * ξ__
    return σ̲̲, α̲̲, κ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾
    # return triu_vec(σ__), triu_vec(α__), κ, triu_vec(ϵ__), triu_vec(ϵₚ__)
    # return nothing
end

function ContinuumMechanicsBase.predict(
            ψ   ::Bammann1990Modeling{T}, # , S},
            test::AbstractBCJMetalTest{T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat} # , S<:SymmetricTensor{2, 3, T}}
    M = ψ.N + 1
    σ̲̲       = zeros(T, 6)   # deviatoric stress
    ϵ̲̲       = zeros(T, 6)   # total strain
    ϵ̲̲⁽ᵖ⁾    = zeros(T, 6)   # plastic strain
    α̲̲       = fill(1e-7, 6) # kinematic hardening
    κ       = 0.0           # isotropic hardening
    ϵ⃗ = []
    σ⃗ = []
    push!(ϵ⃗, ϵ̲̲)
    push!(σ⃗, σ̲̲)
    for i ∈ range(2, M)
        σ̲̲, α̲̲, κ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾ = update(ψ, σ̲̲, α̲̲, κ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, p)
        # update!(ψ, σ__, α__, κ, ϵ__, ϵₚ__, p)
        push!(ϵ⃗, ϵ̲̲)
        push!(σ⃗, σ̲̲)
    end
    return (data=(ϵ=hcat(ϵ⃗...), σ=hcat(σ⃗...)),)
    # σ__     = zeros(T, 6)   # deviatoric stress
    # ϵ__     = zeros(T, 6)   # total strain
    # ϵₚ__    = zeros(T, 6)   # plastic strain
    # α__     = fill(1e-7, 6) # kinematic hardening
    # κ       = 0.0           # isotropic hardening
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