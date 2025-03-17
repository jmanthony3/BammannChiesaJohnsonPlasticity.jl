export DK
export map, predict

using ContinuumMechanicsBase
using ComponentArrays, StructArrays
# using Tensors # uncomment when we can work with Tensors.jl

"""
Structure for viscoplasticity model with loading conditions and material properties.
Here, uses the effective strain rate based on applied strain rate and loading direction.
"""
struct DK{T<:AbstractFloat} <: AbstractBCJMetalModel
# struct DK{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: AbstractBCJMetalModel
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

Use loading conditions and material properties to construct viscoplasticity model which assumes a Poisson's Ratio of 0.5.
Here, `μ` is the shear modulus.
"""
function DK(conditions::BCJMetalStrainControl, μ::AbstractFloat)
    θ       = conditions.θ
    ϵ̇       = conditions.ϵ̇
    ϵₙ      = conditions.ϵₙ
    N       = conditions.N
    loaddir = conditions.loaddir
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
    return DK{T}(θ, ϵ̇_eff, ϵₙ, N, Δϵ, Δt, μ)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function update(ψ::DK, σ̲̲, α̲̲, κ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, (;
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
    Y  *= (C₁₉ < 0.0) ? (1.0) : (0.5 * ( 1.0 + tanh(max(0., C₁₉ * ( C₂₀ - θ )))))


    # trial guesses
    α_mag       = norm_symvec(α̲̲)
    # α_mag       = norm(α̲̲)
    α_mag      *= sqrt23
    # trial guesses: ISVs (from recovery) and stress
    recovery    = Δt * (r_d * ϵ̇_eff + r_s) * α_mag      # recovery for alpha (kinematic hardening)
    Recovery    = Δt * (R_d * ϵ̇_eff + R_s) * κ          # recovery for kappa (isotropic hardening)
    α̲̲⁽ᵗʳ⁾       = α̲̲   .* (1 - recovery)
    # αₜᵣ__       = α̲̲   * (1 - recovery)
    κ⁽ᵗʳ⁾       = κ     * (1 - Recovery)
    σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ * Δϵ)                         # trial stress
    ξ̲̲⁽ᵗʳ⁾       = σ̲̲⁽ᵗʳ⁾ - (2.0 / 3.0) * α̲̲⁽ᵗʳ⁾           # trial overstress original
    # ξ̲̲⁽ᵗʳ⁾       = σ̲̲⁽ᵗʳ⁾ - α̲̲⁽ᵗʳ⁾                       # trial overstress original
    # ξ̲̲⁽ᵗʳ⁾      .= σ̲̲⁽ᵗʳ⁾ - sqrt23 .* α̲̲⁽ᵗʳ⁾             # trial overstress FIT
    ξ_mag       = norm_symvec(ξ̲̲⁽ᵗʳ⁾)
    # ξ_mag       = norm(ξ__)


    # yield criterion
    flow_rule = ξ_mag - sqrt23 * (κ⁽ᵗʳ⁾ + β)            # same as vumat20
    if flow_rule <= 0.0     # elastic
        # trial guesses are correct
        σ̲̲           = @. σ̲̲⁽ᵗʳ⁾
        α̲̲           = @. α̲̲⁽ᵗʳ⁾
        κ           = κ⁽ᵗʳ⁾
        ϵ̲̲          += @. Δϵ
        # state.ϵ_dot_plastic__    .= 0.
    else                    # plastic
        # Radial Return
        Δγ          = flow_rule / (2μ + 2(h + H) / 3)     # original
        n̂           = ξ̲̲⁽ᵗʳ⁾ ./ ξ_mag
        # n̂           = ξ__ / ξ_mag
        σ̲̲_prev      = σ̲̲
        σ̲̲           = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) .* n̂
        α̲̲           = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) .* n̂
        # σ̲̲           = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) * n̂
        # α̲̲           = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) * n̂
        κ           = κ⁽ᵗʳ⁾   + (H * sqrt23 * Δγ)  # original
        ϵ̲̲⁽ᵖ⁾       += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) ./ 2μ))
        # ϵ̲̲⁽ᵖ⁾       += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) / 2μ))
        ϵ̲̲          += @. Δϵ
    end
    # ϵ_dot_plastic__ .= (f * sinh(V \ (ξ_mag - κ - Y)) / ξ_mag) * ξ__
    return σ̲̲, α̲̲, κ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾
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
    σ̲̲       = zeros(T, 6)   # deviatoric stress
    ϵ̲̲⁽ᵖ⁾    = zeros(T, 6)   # plastic strain
    ϵ̲̲       = zeros(T, 6)   # total strain
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
    # ϵₚ__    = zeros(T, 6)   # plastic strain
    # ϵ__     = zeros(T, 6)   # total strain
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
    #     α__     = fill(1e-7, s) # kinematic hardening
    #     κ       = 0.            # isotropic hardening
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