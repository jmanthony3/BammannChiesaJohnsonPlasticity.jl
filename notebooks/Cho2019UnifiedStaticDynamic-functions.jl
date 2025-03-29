# module B93F

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays, StructArrays
# using Tensors # uncomment when we can work with Tensors.jl
using DocStringExtensions

ContinuumMechanicsBase.I₁(tensor::Vector{<:Real}) = sum(x[[1, 4, 6]])
ContinuumMechanicsBase.I₂(tensor::Vector{<:Real}) = 2.0 \ (  ( I₁(tensor) ^ 2.0 )  -  ( I₁(tensor .^ 2.0) )  )
ContinuumMechanicsBase.I₃(tensor::Vector{<:Real}) = det([
    tensor[1] tensor[2] tensor[3];
    tensor[2] tensor[4] tensor[5];
    tensor[3] tensor[5] tensor[6]
])

"""
Structure for viscoplasticity model with loading conditions and material properties.
Here, uses the effective strain rate based on applied strain rate and loading direction.
"""
struct Cho2019Unified{T<:AbstractFloat} <: BammannChiesaJohnsonPlasticity.AbstractBCJMetalModel
# struct Bammann1993Failure{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: AbstractBCJMetalModel
    θ       ::T         # applied temperature
    P       ::T         # pressure
    ϵ̇_eff   ::T         # strain rate (effective)
    ϵₙ      ::T         # final strain
    N       ::Integer   # number of strain increments
    Δϵ      ::Vector{T} # S         # total strain tensor step
    Δt      ::T         # time step
end

"""
    $(SIGNATURES)

Outer constructor for loading conditions and material properties which assumes a Poisson's ratio of 0.5.
Here, `μ` is the shear modulus.
"""
function Cho2019Unified(Ω::BammannChiesaJohnsonPlasticity.BCJMetalStrainControl, P::AbstractFloat)
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
        half_sqrt_three = 0.5 * √(3.0)
        ϵₙ *= half_sqrt_three
        Δϵ .= [0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0]
        # Δϵ  = S([0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0])
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[2] / (ϵ̇ * half_sqrt_three)         # timestep
        # Δt  = Δϵ[1, 2] / ϵ_dot      # timestep
        ϵ̇
    end
    return Cho2019Unified{T}(θ, P, ϵ̇_eff, ϵₙ, N, Δϵ, Δt)
end

"""
Using the equations and constants from [Bammann et. al. (1993)](@cite bammannFailureDuctileMaterials1993), this kernel function maps the current material state and ISVs onto the next configuration.
Though not explicitly listed in paper, temperature equations `h = C₁₅ * exp(-C₁₆ / θ)` and `H = C₁₇ * exp(-C₁₈ / θ)` are included (and their constants renumbered) from (c. f. [Horstemeyer (1994)](@cite horstemeyerPredictingFormingLimit1994)).
Important: `ϕ` is included in the list of arguments, but is presently, internally set to zero.
This is a limitation of the point simulator causing infinite stress triaxiality, χ.
"""
function update(ψ::Cho2019Unified, σ̲̲, α̲̲, κ, ϕ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, (;
            # BCJ-plasticity
            ## yield surface
            C₁,     C₂,     # V
            C₃,     C₄,     # Y
            C₅,     C₆,     # f
            ## pressure-dependent yield surface
            tanB, Pk1, Pk2, Pk3, B2,
            ## kinematic hardening
            C₇,     C₈,     # r_d
            C₉,     C₁₀,    # h
            C₁₁,    C₁₂,    # r_s
            ## isotropic hardening
            # C₁₃,    C₁₄,    # R_d
            # C₁₅,    C₁₆,    # H
            # C₁₇,    C₁₈,    # R_s
            RD0, RDE, RDV, H0, RS0, RSE, RSV, NK,
            ## torsion, tension/compression
            ca, cb,
            ## dynamic recrystallization
            Cx1, Cx2, Cdp,
            Cx3, Cx4, Csp,
            Cx5, Cxa, Cxb, Cxc,
            ## static RX (grain growth)
            sxn, sk0, sxE, sxV, R,
            ## grain size
            d0, Cg1, Cg2, Cg3, z,
            ## damage
            dd, ff, Kic, aa, bb, cc,
            Cnuc, Tnuc, rr0, nn, Tgrw,
            ## irradiation hardening
            kr1, krt, kr2, kr3, kp1, kpt, kp2
        ))
    # get fields from model
        θ       = ψ.θ
        pres    = ψ.P
        ϵ̇_eff   = ψ.ϵ̇_eff
        μ       = ψ.μ
        ϵₙ      = ψ.ϵₙ
        N       = ψ.N
        Δϵ      = ψ.Δϵ
        Δt      = ψ.Δt
        M       = N + 1
        T       = typeof(float(θ))
    # calculation constants/functions
        sqrt_twothirds = √(2.0 / 3.0)
        sqrt_threehalves = √(3.0 / 2.0)
        p(σ)= sum(σ[[1, 4, 6]]) / 3.0
    # irradiation before damage
        Tirr = pres
        M0, Si, damirr = 0.0, 0.0, 1.0
        if Tirr != 0.0
            kr = kr1 * exp( krt / Tirr )
            Si = (kr * flu) ^ (1.0 / kr2)
            M0 = Kr3 * Si
            damirr = exp(   ( kp1 * exp( kpt / Tirr ) * flu ) ^ ( 1.0 / kp2 )   )
        end
    # pressure-temperature dependent reference density
        #--- Olivine paramters
            ttop    = 300.0
            K0      = 129.0     * 1e3
            dK0dT   = -1e-3     * 1e3
            dKdP    = 4.47
            rho0    = 3345.0
            alp     = 3.5       * 1e-5
            G0      = 79.0      * 1e3
            dG0dT   = -0.014    * 1e3
            dG0dP   = 1.60
            ddG0ddP = -0.04     / 1e3
        #--- Cpx (diopside) paramters (Li & Neuville, 2010)
            # ttop    = 300.0
            # K0      = 113.0     * 1e3
            # dK0dT   = -0.012    * 1e3
            # dKdP    = 4.7
            # rho0    = 3270.0
            # alp     = 3.2       * 1e-5
            # G0      = 40.0      * 1e3 # 73.0 GPa (for real diopside)
            # dG0dT   = -0.011    * 1e3
            # dG0dP   = 1.50
            # ddG0ddP = -0.04     / 1e3
        #--- Py (Pyrope) paramters (Hu et al.,2106)
            # ttop    = 300.0
            # K0      = 169.0     * 1e3
            # dK0dT   = -0.0204   * 1e3
            # dKdP    = 4.31
            # rho0    = 3270.0
            # alp     = 2.724     * 1e-5
            # G0      = 90.0      * 1e3
            # dG0dT   = -0.0126   * 1e3
            # dG0dP   = 1.71
            # ddG0ddP = -0.0415   / 1e3
        #--- B-M: density calculation at given pressure and temperature
            KT0   = K0      + (            dK0dT * ( temp - ttop )      )
            RT0   = rho0    * (   1.0 - (    alp * ( temp - ttop )  )   )
            RRT0  = 1.0
            itmax = 10
            convg = 1e-12
            Niter = 0
            # Newton iterations begin
            for k in range(0,itmax)
                RRT073 = RRT0 ^ (7.0/3.0)
                RRT053 = RRT0 ^ (5.0/3.0)
                RRT023 = RRT0 ^ (2.0/3.0)
                FF = pres - (#=[=#
                    1.5 * KT0 * (RRT073 - RRT053) * (
                            1.0 + (  (3.0 / 4.0) * (dKdP - 4.0) * (RRT023 - 1.0)  )
                        )   #=]=#)
                # Define Derivative of F, dF
                dF1 = (18.0 / 24.0) * (dKdP - 4.0) * KT0 * ( RRT0 ^ (-1.0 / 3.0) ) * (RRT073 - RRT053)
                dF2 = 1.5 * KT0
                dF2 = dF2 * ((7.0 / 3.0) * ( RRT0 ^ (4.0 / 3.0) ) - (5.0 / 3.0) * (RRT0 ^ (2.0 / 3.0)))
                dF2 = dF2 * ((3.0 / 4.0) * (dKdP - 4.0) * (RRT023 - 1.0) + 1.0)
                dF  = -(dF1 + dF2)
                # Update Solution
                RRT0 -= FF/dF
                # Convergence Check
                err = abs(dRRT0)
                err <= convg ? break : Niter += 1
                Niter >= (itmax - 1) ? println("BM convergence issue! ", err) : nothing
            end
            # pressure-temperature dependent reference density
            ρ = RRT0 * RT0
    # shear modulus
        #--- 3rd-order Finite Strain (Birch-Murnaghan EOS)
        KT0 = K0   + (           dK0dT * (θ - ttop)      )
        RT0 = rho0 * (   1.0 - (   alp * (θ - ttop)  )   )
        GT0 = G0   + (           dG0dT * (θ - ttop)      )
        b1  = (3KT0 * dG0dP) - 5GT0
        b2  = 9.0(   (KT0 ^ 2.0) * (  ddG0ddP + ( (1.0 / KT0) * (dKdP - 4.0) * dG0dP )  )
                    + 35.0GT0 / 9.0)
        f   = 0.5(   (  (ρ / RT0) ^ (2.0/3.0)  ) - 1.0   )
        μ   = max(0.01, ((1.0 + 2.0f) ^ 2.5) * (GT0 + b1*f + 0.5*b2*(f ^ 2.0)))
        ν   = 0.3
        G   = (2.0/3.0)*μ*(1.0 + ν)/(1.0 - 2ν)
        if imat == 1        # OFHC Cu (irradation-ISV model)
           μ    = 5.47e4    - (34.1 * θ)
           G  = 70000.0
        elseif imat == 2    # T91 ferritic steel (Barrett et al., 2018)
           μ    = 1.01e5    - (65.0 * θ)
           G  = 170000.0
        elseif imat == 3    # Ti6Al4V (Hukuhara&Sanpei,1993)
           μ    = 4.5e4     - (20.0 * θ)
           G  = 85000.0
        end
    # deviatoric strain and effective strain rate
        davg = p(Δϵ)
        ddd = sqrt_twothirds * norm_symvec(ϵ̇_eff) / Δt
    # trial damage
        dam1 = 1.0 - ϕ
        dam2 = 1.0 - min(1.0, ϕ̇ * Δt / dam1)
    # hydrostatic pressure
        if pres > 0.0
            P = pres
        else
            P_H = (p(σ̲̲) * dam2) + (3.0G * davg * dam1)
        end
    # deviatoric stress and invariants
        di1 = I₁(σ̲̲)
        dj2 = I₂(σ̲̲)
        dj3 = I₃(σ̲̲)
    # temperature dependent constants
        V   = C₁    * exp( -C₂ / θ )
        Y   = C₃    * exp(  C₄ / θ )
        f   = C₅    * exp( -C₆ / θ )
        if dj2 == 0.0
            djr = 1.0 - (ca * (4.0 / 27.0))
            djh = 1.0 + (cb * (4.0 / 27.0))
        else
            djr = 1.0     -     (#=[=#
                                    ca * (   (  4.0 / 27.0  )  -  (  (dj3 ^ 2.0) / (dj2 ^ 3.0)  )   )
                #=]=#) - (#=[=#     cb * (   dj3   /   (  dj2 ^ 1.5 )   )     #=]=#)
            djh = 1.0     +     (#=[=#
                                    ca * (   (  4.0 / 27.0  )  -  (  (dj3 ^ 2.0) / (dj2 ^ 3.0)  )   )
                #=]=#) + (#=[=#     cb * (   dj3   /   (  dj2 ^ 1.5 )   )     #=]=#)
        end
        rd  =           C₇    * exp( -C₈  / θ )   * djr
        h   = max(0.0,  C₉    * μ        * djh)
        rs  =           C₁₁   * exp( -C₁₂ / θ )
        Rd  =           C₁₃   * exp( -C₁₄ / θ )   * djr
        H   = max(0.0,  C₁₅   * μ        * djh)
        Rs  =           C₁₇   * exp(   -( C₁₈ + (1e6 * P * C₂₆) )  /  (R * θ)   )
        Rdc = Rd * (ddd ^ -0.0)
    # C_xd = cₓ₁  * exp(   -( Cx2 + (P * c_dp))  /  θ   )
    # C_xs = cₓ₃  * exp(   -( Cx4 + (P * c_sp))  /  θ   )
    # C_h  = Cx5  * KAlMut * dt
    # β   = ( Y * (1 - ϕ) ) + (  V * (1 - ϕ) * asinh( ϵ̇_eff / f )  )
    # (^) will need this eventually | (v) this is for point simulator
    β   = Y + (V * asinh( ϵ̇_eff / f ))
    # @show V, Y, f, ϕ, β
    # @show r_d, r_s, R_d, R_s, h, H
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
    α_mag      *= sqrt_threehalves
    # α_mag       = norm(α̲̲)
    # trial guesses: ISVs (from recovery) and stress
    recovery    = Δt * (rd * ( sqrt_twothirds * ϵ̇_eff ) + rs) * α_mag     # recovery for alpha (kinematic hardening)
    Recovery    = Δt * (Rd * ( sqrt_twothirds * ϵ̇_eff ) + Rs) * κ         # recovery for kappa (isotropic hardening)
    # σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ .* (1 - ϕ) .* Δϵ) .- ( (ϕ̇ * Δt) / (1 - ϕ) )  # deviatoric stress (trial)
    # (^) will need this eventually | (v) this is for point simulator
    σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ .* Δϵ)                                    # deviatoric stress (trial)
    σ̲̲′⁽ᵗʳ⁾      = σ̲̲⁽ᵗʳ⁾ - (p(σ̲̲⁽ᵗʳ⁾) .* [1, 0, 0, 1, 0, 1])
    α̲̲⁽ᵗʳ⁾       = α̲̲    .* (1 - recovery)
    # α̲̲⁽ᵗʳ⁾       = α̲̲     * (1 - recovery)
    κ⁽ᵗʳ⁾       = κ     * (1 - Recovery)
    ξ̲̲⁽ᵗʳ⁾       = σ̲̲′⁽ᵗʳ⁾ - α̲̲⁽ᵗʳ⁾                      # over-stress (trial)
    ξ_mag       = norm_symvec(ξ̲̲⁽ᵗʳ⁾)
    # ξ_mag       = norm(ξ__)


    # yield criterion
    F = ξ_mag - (κ⁽ᵗʳ⁾ + β) # * (1 - ϕ)
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
        κ       =    κ⁽ᵗʳ⁾ + ( H * Δγ * sqrt_twothirds)
        # χ       = sinh(#=[=#  ( 2(2m̄ - 1) * p(σ̲̲) ) / ( (2m̄ + 1) * vonMises(σ̲̲) )  #=]=#)
        # ϕ       =    1 - (#={=#
        #         1 + (#=[=#   (1 - ϕ) ^ (1 + m̄) - 1   #=]=#) * exp(#=[=#
        #             (#=d̄: =# sqrt_twothirds * ϵ̇_eff) * χ * (1 + m̄) * Δt   #=]=#)
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
            ψ   ::Cho2019Unified{T}, # , S},
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
ContinuumMechanicsBase.parameters(::Cho2019Unified) = (
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

nothing
# end # end of module
