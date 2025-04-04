# module B93F

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays, StructArrays
using LinearAlgebra
# using Tensors # uncomment when we can work with Tensors.jl
using DocStringExtensions

ContinuumMechanicsBase.I₁(x::Vector{<:Real}) = sum(x[[1, 4, 6]])
ContinuumMechanicsBase.I₂(x::Vector{<:Real}) = 2.0 \ (  ( I₁(x) ^ 2.0 )  -  ( I₁(x .^ 2.0) )  )
ContinuumMechanicsBase.I₃(x::Vector{<:Real}) = det([
    x[1] x[2] x[3];
    x[2] x[4] x[5];
    x[3] x[5] x[6]
])

"Maps a scalar onto the volumetric portion of the flat vector representation of a second-rank tensor."
volumetric(x::AbstractFloat)    = x .* [1, 0, 0, 1, 0, 1]

"Returns the scalar, hydrostatic portion from the flat vector representation of a second-rank tensor."
hydrostatic(x::Vector{<:Real})  = I₁(x) / 3.0

"Returns the deviatoric of the flat vector representation of a second-rank tensor."
deviatoric(x::Vector{<:Real})   = x - volumetric(hydrostatic(x))

# """
# Structure for viscoplasticity model with loading conditions and material properties.
# Here, uses the effective strain rate based on applied strain rate and loading direction.
# """
# struct Cho2019Unified{T<:AbstractFloat} <: BammannChiesaJohnsonPlasticity.AbstractBCJMetalModel
# # struct Bammann1993Failure{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: AbstractBCJMetalModel
#     θ       ::T         # applied temperature
#     E⁺      ::T
#     V⁺      ::T
#     R       ::T
#     d₀      ::T
#     Kic     ::T
#     𝒹       ::T
#     𝒻       ::T
#     η₀      ::T
#     R₀      ::T
#     P       ::T         # pressure
#     ϵ̇_eff   ::T         # strain rate (effective)
#     ϵₙ      ::T         # final strain
#     N       ::Integer   # number of strain increments
#     Δϵ̲̲      ::Vector{T} # S         # total strain tensor step
#     Δt      ::T         # time step
# end

"""
    $(SIGNATURES)

Outer constructor for loading conditions and material properties which assumes a Poisson's ratio of 0.5.
Here, `μ` is the shear modulus.
"""
function Cho2019Unified(Ω::BammannChiesaJohnsonPlasticity.BCJMetalStrainControl,
        # n   ::T,
        # ω₀  ::T,
        E⁺  ::T,        # activation energy for grain growth
        V⁺  ::T,        # activation volume for grain growth
        R   ::T,        # gas constant
        d₀  ::T,        # initial grain size
        # z   ::T,
        Kic ::T,        # fracture toughness
        𝒹   ::T,        # average size of second phase particles
        𝒻   ::T,        # volume fraction of second phase particles
        η₀  ::T=0.0,    # initial void nucleation density
        R₀  ::T=0.0,    # initial void radius
        P   ::T=0.0) where {T<:AbstractFloat}
    θ       = Ω.θ
    ϵ̇       = Ω.ϵ̇
    ϵₙ      = Ω.ϵₙ
    N       = Ω.N
    loaddir = Ω.loaddir
    M       = N + 1
    # T       = typeof(float(θ))
    Δϵ̲̲      = zeros(T, 6)       # strain increment
    # S       = SymmetricTensor{2, 3, T}
    # Δϵ      = zero(S) # strain increment
    Δt      = (ϵₙ / N) / ϵ̇

    # state evaluation - loading type
    ϵ̇_eff = if loaddir ∈ (:tension, :compression)    # uniaxial tension/compression
        δϵ  = ϵₙ / N
        Δϵ̲̲ .= [δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ]
        # Δϵ  = S([δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ])
        Δt  = δϵ / ϵ̇            # timestep
        ϵ̇
    elseif loaddir == :torsion                                 # torsion
        # convert equivalent strain to true shear strain
        half_sqrt_three = 0.5√(3.0)
        ϵₙ *= half_sqrt_three
        Δϵ̲̲ .= [0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0]
        # Δϵ  = S([0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0])
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ̲̲[2] / (ϵ̇ * half_sqrt_three)         # timestep
        # Δt  = Δϵ[1, 2] / ϵ_dot      # timestep
        ϵ̇
    end
    return Cho2019Unified{T}(θ, #=n,=# #=ω₀,=# E⁺, V⁺, R, d₀, #=z,=# Kic, 𝒹, 𝒻, η₀, R₀, P, ϵ̇_eff, ϵₙ, N, Δϵ̲̲, Δt)
end

"""
Using the equations and constants from [Bammann et. al. (1993)](@cite bammannFailureDuctileMaterials1993), this kernel function maps the current material state and ISVs onto the next configuration.
Though not explicitly listed in paper, temperature equations `h = C₁₅ * exp(-C₁₆ / θ)` and `H = C₁₇ * exp(-C₁₈ / θ)` are included (and their constants renumbered) from (c. f. [Horstemeyer (1994)](@cite horstemeyerPredictingFormingLimit1994)).
Important: `ϕ` is included in the list of arguments, but is presently, internally set to zero.
This is a limitation of the point simulator causing infinite stress triaxiality, χ.
"""
# σ̲̲, α̲̲, κ, κₛ, ϕ, ..., ϕ̇, ..., ϵ̲̲, ϵ̲̲⁽ᵖ⁾, t
# function update(ψ::Cho2019Unified, Sig, Al, K, Ks, Phi, Nuc, Vod, dPhi, X, XR, XH, Xd, Xs, d, TE, PE, VE, Alm, t, (;
function update(ψ::Cho2019Unified, t, σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, κₛ, Si, ϕ, η, damirr, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d, (;
            # BCJ-plasticity
            ## yield surface
            # base, exponent
            C₁,     C₂,             # V
            C₃,     C₄,             # Y
            C₅,     C₆,             # f
            ## pressure-dependent yield surface
            Pₖ₁, Pₖ₂, Pₖ₃,
            ## kinematic hardening
            # base, exponent, pressure
            C₇,     C₈,     C₂₁,    # r_d
            C₉,     C₁₀,    C₂₂,    # h
            C₁₁,    C₁₂,    C₂₃,    # r_s
            ## isotropic hardening
            # base, exponent, pressure
            C₁₃,    C₁₄,    C₂₄,    # R_d
            C₁₅,    C₁₆,    C₂₅,    # H
            C₁₇,    C₁₈,    C₂₆,    # R_s
                            NK,     # * [20250402T1521] (JMA3): I think this is the modifier for finding the k-root
                                    # *                         (see Eq. 4.22 in HEC dissertation)
                                    # *                         (c. f. `optimize.py` that NK=2.0 by default)
            ## torsion, tension/compression
            ca, cb,
            ## dynamic recrystallization
            Cx1, Cx2, Cdp,
            Cx3, Cx4, Csp,
            Cx5, Cxa, Cxb, Cxc,
            ## static RX (grain growth)
            n, ω₀, # E⁺, V⁺, R,
            ## grain size
            # d₀, Cg1, Cg2, Cg3, z,
            Cg1, Cg2, Cg3, z,
            ## damage
            ### nucleation
            # 𝒹, 𝒻, Kic, a, b, c,
            a, b, c,
            # Cnuc, Tnuc, R₀, nn, Tgrw,
            pCnuc, Tnuc, nn, Tgrw,
            ## irradiation hardening
            kr1, krt, kr2, kr3, kp1, kpt, kp2
        ); imat=0, iYS=0, tanβ₀=0.0, iREXmethod=0, iGSmethod=0)
    # get fields from model
        θ       = ψ.θ
        # n       = ψ.n
        # ω₀      = ψ.ω₀
        E⁺      = ψ.E⁺
        V⁺      = ψ.V⁺
        R       = ψ.R
        # z       = ψ.z
        # d₀      = ψ.d₀
        # η₀      = ψ.η₀
        Kic     = ψ.Kic
        𝒹       = ψ.𝒹
        𝒻       = ψ.𝒻
        R₀      = ψ.R₀
        pres    = ψ.P
        P, P_H  = 0.0, 0.0
        ϵ̇_eff   = ψ.ϵ̇_eff
        # μ       = ψ.μ
        ϵₙ      = ψ.ϵₙ
        N       = ψ.N
        # DE0     = ψ.Δϵ
        Δϵ      = ψ.Δϵ̲̲
        ϵ̲̲′      = deviatoric(ϵ̲̲)
        ϵ̲̲′_mag  = norm_symvec(ϵ̲̲′)
        # ϵ̲̲′⁽ᵖ⁾   = deviatoric(ϵ̲̲⁽ᵖ⁾)
        ϵ̲̲⁽ᴴ⁾    = hydrostatic(ϵ̲̲)
        # dt      = ψ.Δt
        Δt      = ψ.Δt
        t      += Δt # ! update state variable
        M       = N + 1
        T       = typeof(float(θ))
    # calculation constants/functions
        sqrt_twothirds = √(2.0/3.0)
        sqrt_threehalves = √(3.0/2.0)
    # * [20250402T1345] (JMA3): Moved to `predict`
    # # irradiation before damage
    #     Tirr = pres
    #     M0, Si, damirr = 0.0, 0.0, 1.0
    #     if Tirr != 0.0
    #         kr = kr1 * exp(krt/Tirr)
    #         Si = (kr*flu) ^ (1.0/kr2)
    #         M0 = Kr3 * Si
    #         damirr = exp(  ( kp1 * exp(kpt/Tirr) * flu )  ^  ( 1.0 / kp2 )  )
    #     end
    # pressure-temperature dependent reference density
        #--- Olivine paramters
            ttop    = 300.0
            κ₀      = 129.0     * 1e3
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
            KT0   = κ₀      + (            dK0dT * (θ-ttop)    )
            RT0   = rho0    * (  1.0  -  (   alp * (θ-ttop) )  )
            RRT0  = 1.0
            itmax = 10
            convg = 1e-12
            Niter = 0
            # Newton iterations begin
            for k in range(0, itmax)
                RRT073 = RRT0 ^ (7.0/3.0)
                RRT053 = RRT0 ^ (5.0/3.0)
                RRT023 = RRT0 ^ (2.0/3.0)
                FF = pres    -    (#={=#   1.5KT0   *   (  RRT073  -  RRT053  )   *   (#=[=#
                    1.0  +  ( (3.0/4.0) * (dKdP-4.0) * (RRT023-1.0) )  #=]=#)   #=}=#)
                # Define Derivative of F, dF
                dF1 = ( 18.0 / 24.0 )  *  ( dKdP - 4.0 )  *  KT0  *  ( RRT0 ^ (-1.0/3.0) )  *  ( RRT073 - RRT053 )
                dF2 = 1.5KT0
                dF2 = dF2 * (  ( 7.0 / 3.0 )  *  ( RRT0 ^ (4.0/3.0) )  -  ( 5.0 / 3.0 )  *  ( RRT0 ^ (2.0/3.0) )  )
                dF2 = dF2 * (  ( 3.0 / 4.0 )  *  ( dKdP - 4.0 )  *  ( RRT023 - 1.0 )  +  1.0)
                dF  = -(dF1+dF2)
                # find Corrector
                dRRT0 = -FF/dF
                # update Solution
                RRT0 -= dRRT0
                # convergence Check
                err = abs(dRRT0)
                err <= convg ? break : Niter += 1
                Niter >= (itmax - 1) ? println("BM convergence issue! ", err) : nothing
            end
            # pressure-temperature dependent reference density
            ρ = RRT0 * RT0
    # shear modulus
        #--- 3rd-order Finite Strain (Birch-Murnaghan EOS)
        KT0 = κ₀     +   (          dK0dT * (θ-ttop)    )
        RT0 = rho0   *   (  1.0 - (   alp * (θ-ttop) )  )
        GT0 = G0     +   (          dG0dT * (θ-ttop)    )
        b1  = (3KT0*dG0dP) - 5GT0
        b2  = 9.0(   (  KT0  ^  2.0  )   *   (  ddG0ddP  +  (
                (1.0/KT0) * (dKdP-4.0) * dG0dP )  )   +    (  35.0GT0  /  9.0  )   )
        F   = 0.5(  ( (ρ/RT0) ^ (2.0/3.0) )  -  1.0  )
        μ   = max(  0.01,  ( (1.0+2.0F) ^ 2.5 )  *  ( GT0 + (b1*F) + 0.5b2 * (F^2.0) )  )
        ν   = 0.3
        K   = (2.0/3.0) * μ * (1.0+ν) / (1.0-2ν)
        if     imat == 1    # OFHC Cu (irradation-ISV model)
           μ = 5.47e4    - (34.1*θ)
           K = 70000.0
        elseif imat == 2    # T91 ferritic steel (Barrett et al., 2018)
           μ = 1.01e5    - (65.0*θ)
           K = 170000.0
        elseif imat == 3    # Ti6Al4V (Hukuhara&Sanpei,1993)
           μ = 4.5e4     - (20.0*θ)
           K = 85000.0
        end
    # deviatoric strain and effective strain rate
        # davg    = hydrostatic(Δϵ)
        Δϵ̲̲⁽ᴴ⁾   = hydrostatic(Δϵ)
        # ϵ̲̲′⁽ᵗʳ⁾
        # DE      = deviatoric(Δϵ)
        Δϵ̲̲′     = deviatoric(Δϵ)
        # ddd     = sqrt_twothirds * norm_symvec(Δϵ̲̲′) / Δt
        ϵ̲̲̇′_mag  = sqrt_twothirds * norm_symvec(Δϵ̲̲′) / Δt
    # trial damage
        # dam1 = 1.0 - ϕ
        # dam2 = 1.0 - min(1.0, ϕ̇*Δt/dam1)
        ϕ₁⁽ᵗʳ⁾ = 1.0 - ϕ
        ϕ₂⁽ᵗʳ⁾ = 1.0 - min(1.0, ϕ̇*Δt/ϕ₁⁽ᵗʳ⁾)
    # hydrostatic pressure
        if pres > 0.0
            P = pres
        else
            P_H = (hydrostatic(σ̲̲)*ϕ₂⁽ᵗʳ⁾) + (3.0K*Δϵ̲̲⁽ᴴ⁾*ϕ₁⁽ᵗʳ⁾)
        end
    # deviatoric stress and invariants
        # σ̲̲′
        # S  = deviatoric(Sig)
        σ̲̲′  = deviatoric(σ̲̲)
        ds = σ̲̲′
        di1 = σ̲̲[1] + σ̲̲[4] + σ̲̲[6]
        # di1  = 3.0 * P_H
        # [
        #     26.643521871951265,       # [1]: σ_11
        #     1.2074287379946301e-8,    # [2]: σ_12
        #     1.2074287379946301e-8,    # [3]: σ_13
        #     -13.321760917864202,      # [4]: σ_22
        #     1.2074287379946301e-8,    # [5]: σ_23
        #     -13.321760917864202,      # [6]: σ_33
        # ]
        # [
        #     26.643521871951265,       # [0]: σ_11
        #     -13.321760917864202,      # [1]: σ_22
        #     -13.321760917864202,      # [2]: σ_33
        #     1.2074287379946301e-8,    # [3]: σ_23
        #     1.2074287379946301e-8,    # [4]: σ_13
        #     1.2074287379946301e-8,    # [5]: σ_12
        # ]
        # [26.643521871951265,    -13.321760917864202,   -13.321760917864202,   1.2074287379946301e-8, 1.2074287379946301e-8, 1.2074287379946301e-8,]
        dj2 = 0.5((sum(ds[[1, 4, 6]] .^ 2.0))
            + 2.0(sum(ds[[2, 3, 5]] .^ 2.0)))
        dj3 =   (ds[1]*(ds[4]*ds[6]-ds[3]*ds[3])
                - ds[2]*(ds[2]*ds[6]-ds[3]*ds[5])
                + ds[5]*(ds[2]*ds[3]-ds[4]*ds[5]))
        # di1 = I₁(σ̲̲)
        # dj2 = I₂(ds)
        # dj3 = I₃(ds)
    # temperature dependent constants
        V   = C₁ * exp(-C₂/θ)
        Y   = C₃ * exp( C₄/θ)
        f   = C₅ * exp(-C₆/θ)
        if dj2 == 0.0
            djr = 1.0 - ( ca * (4.0/27.0) )
            djh = 1.0 + ( cb * (4.0/27.0) )
        else
            djr = 1.0    -    (#=[=#
                                    ca   *   (  ( 4.0 / 27.0 )  -  ( (dj3^2.0) / (dj2^3.0) )  )
                #=]=#) - (#=[=#     cb   *   (  dj3  /  ( dj2 ^ 1.5 )  )   #=]=#)
            djh = 1.0    +    (#=[=#
                                    ca   *   (  ( 4.0 / 27.0 )  -  ( (dj3^2.0) / (dj2^3.0) )  )
                #=]=#) + (#=[=#     cb   *   (  dj3  /  ( dj2 ^ 1.5 )  )   #=]=#)
        end
        rd  =           C₇    * exp(    -C₈                  /        θ    )   *   djr
        h   = max(0.0,  C₉    * μ                                              *   djh)
        rs  =           C₁₁   * exp(    -C₁₂                 /        θ    )
        Rd  =           C₁₃   * exp(    -C₁₄                 /        θ    )   *   djr
        H   = max(0.0,  C₁₅   * μ                                              *   djh)
        Rs  =           C₁₇   * exp(  -( C₁₈ + (1e6P*C₂₆) )  /  ( R * θ )  )
        Rdc = Rd * (ϵ̲̲̇′_mag^-0.0)
    # yield surface parameters
        #iYS: 0-Pressure insensitive (Mises);
        #     1-Pressure sensitive (Shear-Mises);
        #     2-Pressure sensitive (TANH)
        if     iYS == 0
            Yₚ = 0.
        elseif iYS == 1
            tanB = tanβ₀
            Pa = Pₖ₁   *   (  ( 1.0 + exp(-Pₖ₂/θ) )  ^  ( -Pₖ₃ )  )
            Pc = 0.0Pa
            Pd = Pa - Pc
            imode, Yₚ = if P <= Pa
                imode, Ft = (P<Pc)   ?   (1, 0.0)   :   (
                    2,  ( 0.5Pd )  *  ( (P-Pc) ^ 2.0 )  *  tanB   )
                Yₚ = (P*tanB) - Ft
                (   imode,                             Yₚ   )
            else
                (       3,   (  Pa  -  0.5Pd  )   *  tanB   )
            end
        elseif iYS == 2
            #Yp = Pk1*exp(-Pk2*θ)*tanh(B2*P[i])
            #Yp = (Pk1*(1. + exp(-Pk2/θ))^(-Pk3))*tanh(B2*P[i])
            β₁ = max(  1e-10,    ( Pₖ₁ - (Pₖ₂*θ) )  )
            Yₚ = β₁  *  tanh( (Pₖ₃/β₁) * P )            
        else
            error("iYS > 2 which is not supported.")
        end
    # viscous stress
        #Be = V*log((ddd + sqrt(ddd^2 + F^2))/F)
        #... Using sinh^n for strain rate-stress curve's smooth connection
        Be = V      *      log(     f     \     (#={=#
                (   ϵ̲̲̇′_mag   ^   (  1.0  /  NK  )   )    +    sqrt(#=[=#
                    (  ( ϵ̲̲̇′_mag ^ (1.0/NK) )  ^  2.0  )  +  (  f  ^  2.0  )   #=]=#)
            #=}=#)     )
    # previous alpha magnitude
        # α̲̲_mag
        # `\lvboxline`  : ⎸
        # `\rvboxline`  : ⎹
        # `\mid`        : ∣
        # `\Vert`       : ‖
        α̲̲_mag   = norm_symvec(α̲̲)
        α̲̲_mag  *= sqrt_threehalves
    # REX Model
        ## REX calculation: separated DRX and SRX equations
            if     iREXmethod == 0 # Euler Method (explicit)
                KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                dAlpha  = (  h  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*ϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  *  (     κ ^ 2.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ  \  ( dKappa + dAlpha )
                Cxd     = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ϵ̲̲̇′_mag  *  Δt  )
                Cxs     = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  Δt  )
                Ch      = Cx5 * KAlMu1 * Δt
                # pX0     = Cd + Cs
                # dXR     = pX0*(X[i-1]^Cxa)*(1. - X[i-1])^Cxb
                # dXH     = Ch*X[i-1]^Cxc
                # dX      = dXR - dXH
                # new trial
                dXd     = Cxd  *  ( X ^ Cxa )  *  ( (1.0-X) ^ Cxb )
                dXs     = Cxs  *  ( X ^ Cxa )  *  ( (1.0-X) ^ Cxb )
                dXR     = dXd + dXs
                dXH     = Ch  *  ( X ^ Cxc )
                dX      = dXR - dXH
                # ? [20250402T1149] (JMA3): Maybe this (v) should be included? It's not originally...
                # * ========================================================================
                # * [20250402T1151] (JMA3): Maybe this section of reassignment is redundant
                # *                         since these get updated at the end anyway.
                # XR     += dXR # ! update ISV
                XH     += dXH # ! update ISV
                Xd     += dXd # ! update ISV
                Xs     += dXs # ! update ISV
                # * ========================================================================
                # dX      = dXR - dXH
                xx      = X + dX
            elseif iREXmethod == 1 # explicit exponential integration algorithm
                KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                dAlpha  = (  h  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*ϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  *  (      κ ^ 2.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ  \  ( dKappa + dAlpha )
                Cxd      = Cx1   *   exp(  -( Cx2 + (   P*Cdp) )  /        θ  )   *   (  KAlMu  *  ϵ̲̲̇′_mag  *  Δt  )
                Cxs      = Cx3   *   exp(  -( Cx4 + (   P*Csp) )  /        θ  )   *   (  KAlMu          *  Δt  )
                Ch      = Cx5 * KAlMu1 * Δt
                Udt     = ( Cxd + Cxs )  *  ( X ^ Cxa )  *  ( (1.0-X) ^ (Cxb-1.0) )
                Udt     = Udt   +   (  Ch  *  (  X ^ (Cxc-1.0)  )  )
                Vdt     = ( Cxd + Cxs )  *  ( X ^ Cxa )  *  ( (1.0-X) ^ (Cxb-1.0) )
                xx      = (  X   *   exp( -Udt )  )   +   (  Vdt  *  ( (1.0-exp(-Udt)) / Udt )  )
                pX0     = Cxd + Cxs
                dXR     = pX0          *  ( X ^ Cxa )  *  ( (1.0-X) ^  Cxb  )
                dXH     = Ch * (X^Cxc)
            elseif iREXmethod == 2 # RK4-explicit method
                # K = 10.
                KAlMu   = μ  \  ( (κ^1.0) + (α̲̲_mag^1.0) )
                dAlpha  = (  h  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*ϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 3.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  *  (      κ ^ 3.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ  \  ( dKappa + dAlpha )
                Cxd      = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /  ( R * θ )  )   *   (  KAlMu  *  ϵ̲̲̇′_mag  )
                # Cd     = Cx1*exp(-(Cx2 + P[i]*Cdp)/θ)*KAlMu*ddd
                Cxs      = Cx3   *   exp(  -( Cx4 + (   P*Csp) )  /        θ    )   *      KAlMu
                Ch      = Cx5 * KAlMu1
                CC      = 0.5(Cxd + Cxs)
                k₁      = CC   *   (  X ^ Cxa  )   *   (  ( 1.0-X ) ^ Cxb  )
                k₁      = k₁ - 0.5(   Ch   *   (    X                 ^  Cxc  )   )
                k₂      = CC    *    (   (  X  +     ( Δt * k₁ )  )   ^   Cxa   )    *    (   (  1.0  -  ( X +    (Δt*k₁) )  )   ^   Cxb   )
                k₂      = k₂ - 0.5(   Ch   *   (  ( X +    (Δt*k₁) )  ^  Cxc  )   )
                k₃      = CC    *    (   (  X  +     ( Δt * k₂ )  )   ^   Cxa   )    *    (   (  1.0  -  ( X +    (Δt*k₂) )  )   ^   Cxb   )
                k₃      = k₃ - 0.5(   Ch   *   (  ( X +    (Δt*k₂) )  ^  Cxc  )   )
                k₄      = CC    *    (   (  X  +  2.0( Δt * k₃ )  )   ^   Cxa   )    *    (   (  1.0  -  ( X + 2.0(Δt*k₃) )  )   ^   Cxb   )
                k₄      = k₄ - 0.5(   Ch   *   (  ( X + 2.0(Δt*k₃) )  ^  Cxc  )   )
                xx      = X   +   (  1.0  /  3.0  )   *   (  ( k₁ + 2.0(k₂+k₃) + k₄ )  *  Δt  )
                Cxd     *= Δt
                Cxs     *= Δt
                Ch     *= Δt
                pX0     = Cxd + Cxs
                dXd     = Cxd  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXs     = Cxs  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXR     = pX0 *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXH     = Ch * (xx ^ Cxc)
            elseif iREXmethod >= 3 # implicitly solve functions using Newton-Rapson method
                Nitmax = 20
                Ntol   = 1e-6
                xx     = 0.5
                for k in range(0, Nitmax)
                    if     iREXmethod == 3 # Euler Method (implicit)
                        KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                        dAlpha  = (  h  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*ϵ̲̲̇′_mag) + rs )  *   ( α̲̲_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (  H  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  *   (      κ ^ 2.0 )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ  \  ( dKappa + dAlpha )
                        Cxd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ϵ̲̲̇′_mag  *  Δt  )
                        Cxs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  Δt  )
                        Ch      = Cx5 * KAlMu1 * Δt
                        # TODO [20250331T1119] (JMA3): come back to decrement this section instead
                        F       = X  +   (  ( Cxd + Cxs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )  )
                        F       = F  - ( Ch * (xx^Cxc) ) - xx
                        dF      =    - ( Cxd + Cxs )  *  Cxb  *  ( (1.0-xx) ^ (Cxb-1.0) )  *  ( xx ^  Cxa )
                        dF      = dF + ( Cxd + Cxs )  *  Cxa  *  ( (1.0-xx) ^  Cxb      )  *  ( xx ^ (Cxa-1.0))
                        dF      = dF - (      Ch    *  Cxc  *  ( xx ^ (Cxc-1.0) )  )
                        dF      = dF - 1.0
                    elseif iREXmethod == 4 # exponential integration algorithm (asymptotic)
                        KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                        dAlpha  = (  h  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*ϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (  H  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  *  (      κ ^ 2.0 )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ  \  ( dKappa + dAlpha )
                        Cxd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ϵ̲̲̇′_mag  *  Δt  )
                        Cxs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  Δt  )
                        Ch      = Cx5 * KAlMu1 * Δt
                        Udt     = ( Cxd + Cxs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        Udt     = Udt    +    (  Ch  *  ( xx ^ (Cxc-1.0) )  )
                        Vdt     = ( Cxd + Cxs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        F       = (   X * exp(-Udt)   )    +    (   Vdt   *   (  ( 1.0-exp(-Udt) )  /  Udt  )   )    -    xx
                        dUdt    = ( Cxc - 1.0 )  *  Ch  *  ( xx ^ (Cxc-2.0) )
                        dUdt    = dUdt   +   (  ( Cxd + Cxs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dUdt    = dUdt   -   (  ( Cxd + Cxs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        dVdt    =            (  ( Cxd + Cxs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dVdt    = dVdt   -   (  ( Cxd + Cxs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        dF      = -X * dUdt * exp(-Udt)
                        dF      = dF    +    (   (  ( dVdt/Udt )  -  ( (Vdt*dUdt) / (Udt^2.0) )  )   *   (   1.0  -  exp(-Udt)  )   )
                        dF      = dF    +    (   (                      Vdt       /  Udt         )   *   (  dUdt  *  exp(-Udt)  )   )
                        dF     -= 1.0
                    elseif iREXmethod == 5 # exponential integration algorithm (trapezoidal)
                        KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                        dAlpha  = (  h  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*ϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (  H  *  ϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  *  (      κ ^ 2.0 )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ  \  ( dKappa + dAlpha )
                        Cxd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ϵ̲̲̇′_mag  *  Δt  )
                        Cxs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  Δt  )
                        Ch      = Cx5 * KAlMu1 * Δt
                        Udt     = ( Cxd + Cxs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        Udt     = Udt   +   (  Ch  *  ( xx ^ (Cxc-1.0) )  )
                        U0dt    = ( Cxd + Cxs )  *  (  X ^ Cxa )  *  ( (1.0- X) ^ (Cxb-1.0) )
                        U0dt    = U0dt  +   (  Ch  *  (  X ^ (Cxc-1.0) )  )
                        Vdt     = ( Cxd + Cxs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        V0dt    = ( Cxd + Cxs )  *  (  X ^ Cxa )  *  ( (1.0- X) ^ (Cxb-1.0) )
                        F       =                   X  *  exp( -0.5(U0dt+Udt) )
                        F      +=          0.5(  V0dt  *  exp( -0.5(U0dt+Udt) )  +  Vdt  )   -   xx
                        dUdt    = ( Cxc - 1.0 )  *  Ch  *  ( xx ^ (Cxc-2.0) )
                        dUdt    = dUdt   +   (  ( Cxd + Cxs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dUdt    = dUdt   -   (  ( Cxd + Cxs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        dVdt    =            (  ( Cxd + Cxs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dVdt    = dVdt   -   (  ( Cxd + Cxs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        dF      = -0.5dUdt  *  X  *  exp( -0.5(U0dt+Udt) )
                        dF     += 0.5(  -0.5dUdt  *  V0dt  *  exp( -0.5(U0dt+Udt) )  +  dVdt  )
                        dF     -= 1.0
                    end

                    dxx = -F / dF
                    xx  = max(1e-6, min(0.9999999, xx + dxx))

                    pX0 =  Cxd  +  Cxs
                    dXR = pX0  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                    dXH =  Ch  *  ( xx ^ Cxc )

                    if abs(dxx) <= Ntol
                        break
                    end
                    if k >= Nitmax-1
                        println("Newton-Rapson Convergence Issue: k >= Nitmax")
                    end
                end
            else
                error("iREXMethod > 5 which is not supported.")
            end

            # Final solution and estimate volume of DRX and SRX
            xx  = max(1e-6, min(0.9999999, xx))
            dX  = xx - X
            X0  = 1.0  -  ( dX / (1.0-X) )
            X   = xx # ! update ISV
            (dX < 0.0)  ?  (X0 = 1.0)  :  nothing # ∵ REX is irreversible
            XR += dXR # ! update ISV
            XH += dXH # ! update ISV
            #dXd   = dXR*(Cd/pX0)
            #dXs   = dXR*(Cs/pX0)
            Xd += dXd # ! update ISV
            Xs += dXs # ! update ISV
        ## Grain size kinetics (SGG and grain refinement rate)
            # Grain size rate integration method
            # 0-explicit; 1-implicit; 2-analytic; 3-earlier model (IJP,2019)
            # iGSmethod = 0 # [20250402T1523] (JMA3): I commented this out to let the positional argument have precedence
            dim1 = d
            if     iGSmethod == 0 # Forward Euler (explicit)
                # static grain growth rate
                dr      = dim1
                dsgk    =  ω₀   *   exp(  -( E⁺ + (1e6P*V⁺) )  /  ( R * θ )  )
                dsgg    = dsgk   /   (  n  *  ( dr ^ (n-1.0) )  )
                # dynamic grain size reduction rate (new version: EPSL2020)
                dred    = Cg1 * X * ϵ̲̲̇′_mag * (dr^Cg2)
                # total grain size change rate
                d       = dr  +  ( (dsgg-dred) * Δt ) # ! update ISV
                # Z       = ddd*exp((sxE + P[i]*1.e6*sxV)/(R*θ))
                # dss     = (sxk/(Cg3*sxn*0.3))^(1./(sxn-1.+Cg2))*Z^(-(1./(sxn-1.+Cg2)))
            elseif iGSmethod == 1 # Backward Euler: a = 1 (implicit); a = 0.5 (Crank-Nicholson)
                λ       = 1.0
                Nitmax  = 20
                Convg   = 1e-6
                dr      = dim1
                # dsgk    = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(R*θ))
                # time downscaling factor for matching to n=4
                tscl    = t  ^  ( (n/4.0) - 1.0 )
                dsgk    = ω₀   *   exp(  -( E⁺ + (1e6P*V⁺) )  /  ( R * θ )  )   *   tscl
                xx      = dr
                for k in range(0, Nitmax)
                    F   = dr      +      (#={=#     dsgk     *     Δt     /     (#=[=#    n    *    (
                            (  ( 1.0 - λ )  *  ( dr ^ (n-1.0) )  )   +   (  λ  *  ( xx ^ (n-1.0) )  )
                        )    #=]=#)     #=}=#)
                    F  -= Cg1   *   X   *   ϵ̲̲̇′_mag   *   Δt   *   (
                            ( (1.0-λ) * (dr^Cg2) )  +  ( λ * (xx^Cg2) )  )
                    F  -= xx
                    dF  = (  ( dsgk * Δt * λ * (1.0-n) / n )  *  ( xx ^ -n )  )   -   1.0
                    dF -= Cg1    *    X    *    ϵ̲̲̇′_mag    *    Δt    *    (
                            (  Cg2 * λ * ( xx ^ (Cg2-1.0) )  )   )
                    dxx = -F / dF
                    xx += dxx

                    if abs(dxx) <= Convg
                        break
                    end
                    if k >= Nitmax - 1
                        println("N-R Convg Issue for Grain Size: k >= Nitmax", dxx, k)
                    end
                end
                d       = xx # ! update ISV
                prefct  = ( ω₀ * tscl / (Cg1*n*X) )  ^  ( 1.0 / (n-1.0+Cg2) )
                dsss    = prefct     *     (    (#=[=#
                        ϵ̲̲̇′_mag   *   exp(  ( E⁺ + (1e6P*V⁺) ) / ( R * θ )  )
                    #=]=#)    ^    (   -1.0   /   (  n  -  1.0  +  Cg2  )   )    )
            elseif iGSmethod == 2 # analytical solution
                # static grain growth
                dsgk    = ω₀   *   exp(  -( E⁺ + (1e6P*V⁺) )  /  ( R * θ )  )
                # ! update ISV
                # ? [20250401T1206] (JMA3): what is `d0`
                d       = ψ.d₀    +    (   dsgk   *   t   *   (  t  ^  ( (n/4.0) - 1.0 )  )   )    ^    (   1.0   /   n   )
            elseif iGSmethod == 3 # original version of DRX grain size kinetics model
                P1      = 300.0
                P2      = 0.18
                P3      = 2.0
                dr      = dim1
                tscl    = t  ^  ( (n/4.0) - 1.0 )
                dsgk    = ω₀   *   exp(  -( E⁺ + (1e6P*V⁺) )  /  ( R * θ )  )   *   tscl
                dssmax  = ( (dsgk*Δt) + (dr^n) )  ^  ( 1.0 / n )
                # ? [20250331T1347] (JMA3): what even is this if-statement?
                if ϵ̲̲̇′_mag * Δt == 0.0
                    dssr = ( (dsgk*Δt) + (dr^n) )  ^  ( 1.0 / n )
                    dssr = dr
                else
                    dss0 = ϵ̲̲̇′_mag   *   exp(  ( E⁺ + (1e6P*V⁺) )  /  ( R * θ )  )
                    dssr = P1 * (dss0^-P2)
                end
                # ? [20250331T1350] (JMA3): why the addition, subtraction, and increment?
                ddgrw   = (  ( (dsgk*Δt) + (dr^n) )  ^  ( 1.0 / n )  )   -   dr
                dr     += dr + ddgrw
                dss     = min(dssr, dr)
                ddred   = -P3 * X * ϵ̲̲̇′_mag * Δt * dr * (dr-dss)
                # ddred = -Cg3*Xd[i]*dr*(dr - dss)*ddd*dt
                ddred   = max((dss-dr), ddred)
                d       = dr + ddred # ! update ISV
            else
                error("iGSmethod > 3 not supported")
            end
        ## Hall-Petch effect
            idzz = 0
            # ? [20250401T1206] (JMA3): what is `d0`
            dzz1, dzz0 = if idzz == 0
                ( (ψ.d₀/d) ^ z,            1.0 )
            elseif idzz == 1
                (         1.0,   (dim1/d) ^ z )
            elseif idzz == 2
                ( (ψ.d₀/d)     ,   (dim1/d)      ) .^ z
            else
                error("idzz > 2 which is not supported.")
            end
            # d0 = 1. !Turn on if absolute grain size-stress relation is used
            # YT  = YT*dzz1 ! Turn on if grain size dependent yield is used
    # elastic prediction
        twoμ = 2.0μ
        #--- trial deviatoric stress
            # for k in range(0, 3)
            #     Se[k]  = twoμ*DE[k]
            #     Str[k] = (S[k][i-1]*dam2) + (Se[k]*dam1)
            # end
            # for k in range(3, 6)
            #     Se[k]  = twoμ*DE[k]
            #     Str[k] = (S[k][i-1]*dam2) + (Se[k]*dam1)
            # end
            # σ̲̲′⁽ᵗʳ⁾
            # Str = ( σ̲̲′ .* ϕ₂⁽ᵗʳ⁾ )  +  ( twoμ .* Δϵ̲̲′ .* ϕ₁⁽ᵗʳ⁾)
            σ̲̲′⁽ᵗʳ⁾ = ( σ̲̲′ .* ϕ₂⁽ᵗʳ⁾ )  +  ( twoμ .* Δϵ̲̲′ .* ϕ₁⁽ᵗʳ⁾)
        #--- use of Newton Method for DG and Kappa
        iNewton = 0
        #--- irradiation hardening effect
        Hir = (1.0+Si) ^ 2.0
        #--- trial kappa
        if iNewton == 0
            rdrsk   = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) )  *  Δt  *  ( κ ^ (NK-1.0) )  *  dzz1  )
            # Ktr     = κ * X0 * dzz0 / rdrsk
            κ⁽ᵗʳ⁾   = κ * X0 * dzz0 / rdrsk
        end
        #--- trial M in isotropic hardening (output only)
        rdrssk  = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) )  *  Δt  *  ( κₛ ^ (NK-1.0) ) *  dzz1  )
        # Kstr    = κₛ * X0 * dzz0 / rdrssk
        κₛ⁽ᵗʳ⁾  = κₛ * X0 * dzz0 / rdrssk
        # ! update ISV
        # ? [20250401T1206] (JMA3): why are we updating this again?
        d       =           (  ( Rs + (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) )  *         ( κₛ ^ (NK-1.0) )          )
        if iNewton == 1 # Newton iteration (Backward Euler)
            Nitmax  = 20
            Ntol    = 1.e-06
            Rx0     = X0
            xx      = κ
            for k in range(0, Nitmax)
                RSRD    = 1.0  +  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs ) * Δt * ( xx ^ (NK-1.0) )
                F1      = (Rx0*κ/RSRD) - xx
                dF1     = (  -Rx0  *  κ  *  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  *  Δt  *  ( NK - 1.0 )  *  ( xx ^ (NK-2.0) )  /  ( RSRD ^ 2.0 )  )   -   1.0
                dxx     = -F1 / dF1
                xxn     = xx
                xx     += dxx
                if abs(dxx / xxn) <= Ntol
                    break
                end
                if k >= Nitmax - 1
                    println("Ktr: N-R Conv Issue: k >= Nitmax", k)
                end
            end
            κ⁽ᵗʳ⁾   = xx
            rdrsk   = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) )  *  Δt  *  ( xx ^ (NK-1.0) )  )
        end
        #--- trial alpha
        # * [20250331T1452] (JMA3): this is already defined above
        # Al_mag = Al[0][i-1]^2 + Al[1][i-1]^2 + Al[2][i-1]^2 \
        #        +(Al[3][i-1]^2 + Al[4][i-1]^2 + Al[5][i-1]^2)*2.
        # Al_mag = sqrt(Al_mag)*sqrt_threehalves
        rdrsa   = 1.0   +   (  ( rs + (sqrt_twothirds* rd*ϵ̲̲̇′_mag) )  *  Δt  *  α̲̲_mag            *  dzz1  )
        # for k in range(0, 6)
        #     Altr[k] = Al[k][i-1] * X0 * dzz0 / rdrsa
        # end
        # α̲̲⁽ᵗʳ⁾
        # Altr = α̲̲ .* X0 .* dzz0 ./ rdrsa
        α̲̲⁽ᵗʳ⁾ = α̲̲ .* X0 .* dzz0 ./ rdrsa
        #--- Plastic direction tensor N
        # for k in range(0, 6)
        #     Xi[k] = Str[k] - ( (2.0/3.0) * Altr[k] )
        # end
        # Xi_mag2 = Xi[0]^2 + Xi[1]^2 + Xi[2]^2 \
        #         +(Xi[3]^2 + Xi[4]^2 + Xi[5]^2)*2.0
        # Xi_mag = sqrt(Xi_mag2)
        # for k in range(0, 6)
        #     N[k] = Xi[k] / Xi_mag
        # end
        # ξ̲̲⁽ᵗʳ⁾
        # ξ̲̲⁽ᵗʳ⁾_mag
        # n̂
        # Xi = σ̲̲′⁽ᵗʳ⁾  -  ( (2.0/3.0) .* α̲̲⁽ᵗʳ⁾ )
        # Xi_mag = norm_symvec(Xi)
        # N = Xi ./ Xi_mag
        ξ̲̲′⁽ᵗʳ⁾      = σ̲̲′⁽ᵗʳ⁾  -  ( (2.0/3.0) .* α̲̲⁽ᵗʳ⁾ )
        ξ̲̲′⁽ᵗʳ⁾_mag  = norm_symvec(ξ̲̲′⁽ᵗʳ⁾)
        n̂′          = ξ̲̲′⁽ᵗʳ⁾ ./ ξ̲̲′⁽ᵗʳ⁾_mag
    # check plasticity
        ak     = Y + κ⁽ᵗʳ⁾ + Be + Yₚ
        critra = ξ̲̲′⁽ᵗʳ⁾_mag - (sqrt_twothirds*ak*ϕ₁⁽ᵗʳ⁾)
        #if(TEm[i-1] > 0.02): critra = 100.
    # Radial-Return
    if critra <= 0.0 # elastic solution update
        EPflag = 1
        # deviatoric stress update
        # for k in range(0, 6)
        #     S[k][i]   = Str[k]
        #     Sig[k][i] = S[k][i]
        # end
        # σ̲̲ = @. σ̲̲⁽ᵗʳ⁾
        σ̲̲′ .= σ̲̲′⁽ᵗʳ⁾
        # Cauchy stress update
        # for k in range(0, 3)
        #     Sig[k][i] = Sig[k][i] + P_H
        # end
        σ̲̲ .= σ̲̲′⁽ᵗʳ⁾ + volumetric(P_H) # ! update ISV
        # von Mises stress update
        # vM[i]  = S[0][i]^2 + S[1][i]^2 + S[2][i]^2 \
        #         +(S[3][i]^2 + S[4][i]^2 + S[5][i]^2)*2.
        # vM[i]  = sqrt(vM[i])*sqrt_threehalves
        vM = sqrt_threehalves * norm_symvec(σ̲̲′)
        # kinematic hardening & total strain update
        # for k in range(0, 6)
        #     Al[k][i] = Altr[k]
        #     TE[k][i] = TE[k][i-1] + DE[k]
        # end
        # α̲̲ = α̲̲⁽ᵗʳ⁾
        α̲̲ = α̲̲⁽ᵗʳ⁾ # ! update ISV
        # [20250402T1156] (JMA3):   This (v) may not be needed since `ϵ̲̲′` is evaluated
        #                           at the top and should be updated in `predict()`.
        ϵ̲̲′ += Δϵ̲̲′
        # kinematic hardening update
        # [20250402T1156] (JMA3):   This (v) may not be needed since `ϵ̲̲′` is evaluated
        #                           at the top and should be updated in `predict()`.
        # Alm[i] = Al[0][i]^2 + Al[1][i]^2 + Al[2][i]^2 \
        #         +(Al[3][i]^2 + Al[4][i]^2 + Al[5][i]^2)*2.
        # Alm[i] = sqrt(Alm[i])*sqrt_threehalves
        α̲̲_mag = sqrt_threehalves * norm_symvec(α̲̲)
        # isotropic hardening update
        # K[i]   = Ktr
        κ = κ⁽ᵗʳ⁾ # ! update ISV
        # irradiation hardening update
        # Ks[i]  = Kstr
        κₛ = κₛ⁽ᵗʳ⁾ # ! update ISV
        # total equivalent strain update
        # [20250402T1156] (JMA3):   This (v) may not be needed since `ϵ̲̲′` is evaluated
        #                           at the top and should be updated in `predict()`.
        # TEm[i] = TEm[i-1] + (ddd*dt)
        ϵ̲̲′_mag += ϵ̲̲̇′_mag*Δt
        # volumetric strain update
        # [20250402T1156] (JMA3):   This (v) may not be needed since `ϵ̲̲′` is evaluated
        #                           at the top and should be updated in `predict()`.
        # VE[i]  = VE[i-1] + (3.0*davg)
        ϵ̲̲⁽ᴴ⁾ += 3.0Δϵ̲̲⁽ᴴ⁾
        # damage update
        # Phi[i] = Phi[i-1]
        # Nuc[i] = Nuc[i-1]
        # Vod[i] = Vod[i-1]
        # dPhi[i] = 0.0
        ϕ = ϕ # ! update ISV
        η = η # ! update ISV
        νᵥ = νᵥ # ! update ISV
        ϕ̇ = 0.0 # ! update ISV

        # [20250401T1450] (JMA3): I'm not really sure what this is doing; maybe TERRAfit?
        α̲̲ₛₐₜ_mag    = 0.0
        κₛₐₜ        = (  (H * Hir * ϵ̲̲̇′_mag )  /  ( (sqrt_twothirds*Rdc*ϵ̲̲̇′_mag) + Rs )  )   ^   (  1.0  /  NK  )
        ϵ̲̲̲̇′⁽ᵖ⁾_mag   = ϵ̲̲̇′_mag # [20250402T1207] (JMA3): 
        vMₛₐₜ       = Be + Y + Yₚ + α̲̲ₛₐₜ_mag + κₛₐₜ
    else # plastic solution (Radial return starts)
        EPflag = 2
        #--- Plastic strain increment solution
        if     iNewton == 0 # analytical solution for DG
            Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   ak   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
                (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   ϕ₁⁽ᵗʳ⁾   *   (  2.0  /  3.0  )   *   (
                        ( (1.0-X) ^ NK )  *  dzz1  *  ( (h/rdrsa) + (H*Hir/rdrsk) )  )   )    )
        elseif iNewton == 1 # Newton-Rapson for DG and Kappa
            Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   ak   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
                (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   ϕ₁⁽ᵗʳ⁾   *   (  2.0  /  3.0  )   *   (
                        ( (1.0-X) ^ NK )  *  dzz1  *  ( (h/rdrsa) + (H*Hir/rdrsk) )  )   )    )

            Nitmax  = 20
            Ntol    = 1e-6
            xx1     = Δγ
            xx2     = κ
            κ₀      = κ
            Rx      = (1.0-X) ^ NK
            th      = 1.0 # 1-Backward Euler; 0.5-Midpoint; 0-Forward Euler
            for k in range(0, Nitmax)
                thK0thK     = (  ( 1.0 - th )  *  ( κ₀ ^ (NK-1.0) )  )   +   (  th  *  ( xx2 ^ (NK-1.0) )  )
                Rdxx1Rsdt   = 1.0  +  (  ( (sqrt_twothirds*Rdc*xx1) + (Rs*Δt) )  *  thK0thK  )
                F₁      = ξ̲̲′⁽ᵗʳ⁾_mag   -   (  twoμ * xx1  )   -   (
                        sqrt_twothirds  *  ( κ₀ + (Rx*H*Hir*xx1) )  /  Rdxx1Rsdt  )   -   (
                        sqrt_twothirds  *  ( Be + Y + Yₚ )  )
                ∂F₁╱∂x₁ = -twoμ     -     (    sqrt_twothirds    *    (
                    (  Rx  *  H  *  Hir  /  Rdxx1Rsdt  )   -   (
                        ( κ₀ + (Rx*H*Hir*xx1) )  *  sqrt_twothirds  *  Rdc  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )   )    )
                ∂F₁╱∂x₂ = -1.0 * sqrt_twothirds
                F₂      = (  ( κ₀ + (Rx*H*Hir*xx1) )  /  Rdxx1Rsdt  )   -   xx2
                ∂F₂╱∂x₁ = (  Rx  *  H  *  Hir  /  Rdxx1Rsdt  )   -   (
                    ( κ₀ + (Rx*H*Hir*xx1) )  *  sqrt_twothirds  *  Rdc  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )
                ∂F₂╱∂x₂ = -( κ₀ + (Rx*H*Hir*xx1) )  *  (
                    (sqrt_twothirds*Rdc*xx1) + (Rs*Δt) )  *  th  *  ( NK - 1.0 )  *  ( xx2 ^ (NK-2.0) )
                ∂F₂╱∂x₂ = ( ∂F₂╱∂x₂ / (Rdxx1Rsdt^2.0) )  -  1.0

                a₁₁     = ∂F₁╱∂x₁
                a₁₂     = ∂F₁╱∂x₂
                a₂₁     = ∂F₂╱∂x₁
                a₂₂     = ∂F₂╱∂x₂
                dxx2    = ( (-F₂*a₁₁/a₂₁) + F₁ )  /  ( (a₂₂ * (a₁₁/a₂₁)) - a₁₂ )
                dxx1    = ( (-a₁₂*dxx2)   - F₁ )  /  a₁₁

                xx1n    = xx1
                xx2n    = xx2

                xx1    += dxx1
                xx2    += dxx2

                if abs(dxx1/xx1n) <= Ntol && abs(dxx2/xx2n) <= Ntol
                    break
                end
                if k >= Nitmax - 1
                    println("Gamma-K: N-R Conv. Issue: k >= Nitmax")
                end
            end
            Δγ  = xx1
            κ   = xx2 # ! update ISV
        end
        #--- stress solution
            # deviatoric stress update
            # for k in range(0, 6)
            #     S[k][i] = Str[k] - (dam1 * twoμ * DG * N[k])
            # end
            σ̲̲′ .= σ̲̲′⁽ᵗʳ⁾  -  ( (ϕ₁⁽ᵗʳ⁾*twoμ*Δγ) .* n̂′ )
            # Cauchy stress update
            # for k in range(0, 3)
            #     Sig[k][i] = S[k][i] + P_H
            # end
            # for k in range(3, 6)
            #     Sig[k][i] = S[k][i]
            # end
            σ̲̲ .= σ̲̲′ + volumetric(P_H)
            # von Mises stress update
            # vM[i] = S[0][i]^2 + S[1][i]^2 + S[2][i]^2 \
            #         +(S[3][i]^2 + S[4][i]^2 + S[5][i]^2)*2.
            # vM[i] = sqrt(vM[i])*sqrt_threehalves
            vM = sqrt_threehalves * norm_symvec(σ̲̲′)
        #--- total deviatoric strain
        # for k in range(0, 6)
        #     TE[k][i] = TE[k][i-1] + DE[k]
        # end
        ϵ̲̲′ += Δϵ̲̲′
        #--- total plastic strain
        # PE[i] = PE[i-1] + (sqrt_twothirds * DG)
        ϵ̲̲⁽ᵖ⁾ += ( (sqrt_twothirds*Δγ) .* n̂′ ) # ! update ISV
        #--- total volumetric strain
        # VE[i] = VE[i-1] + (3.0 * davg)
        ϵ̲̲⁽ᴴ⁾ += 3.0Δϵ̲̲⁽ᴴ⁾
        #--- alpha solution
        # for k in range(0, 6)
        #     Al[k][i] = Altr[k]   +   (
        #         ( (1.0-X[i]) ^ NK )  *  dzz1  *  dam1  *  h  *  DG  *  N[k]  /  rdrsa  )
        # end
        α̲̲ .= α̲̲⁽ᵗʳ⁾    +    ( # ! update ISV
                (  ( (1.0-X) ^ NK )  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  h  *  Δγ  )   .*   n̂′   ./   rdrsa   )
        #--- kappa solution # ! update ISV
        κ = (   iNewton   !=   0   )    ?    (   xx2   )    :    (#=[=#   κ⁽ᵗʳ⁾   +   (
                ( (1.0-X) ^ NK )  *  sqrt_twothirds  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  H  *  Hir  *  Δγ  /  rdrsk  )   #=]=#)
        #--- irradiation hardening solution in isotropic hardening
        Sir = Si^2.0
        κₛ  = κₛ⁽ᵗʳ⁾   +   ( # ! update ISV
            ( (1.0-X) ^ NK )  *  sqrt_twothirds  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  H  *  Sir  *  Δγ  /  rdrssk  )
        #--- damage
            # for k in range(0, 3)
            #     ds[k] = S[k][i]
            # end
            # for k in range(3, 6)
            #     ds[k] = S[k][i]
            # end
            ds .= σ̲̲′
            di1 = σ̲̲[1] + σ̲̲[4] + σ̲̲[6]
            # di1  = 3.0 * P_H
            # [
            #     26.643521871951265,       # [1]: σ_11
            #     1.2074287379946301e-8,    # [2]: σ_12
            #     1.2074287379946301e-8,    # [3]: σ_13
            #     -13.321760917864202,      # [4]: σ_22
            #     1.2074287379946301e-8,    # [5]: σ_23
            #     -13.321760917864202,      # [6]: σ_33
            # ]
            # [
            #     26.643521871951265,       # [0]: σ_11
            #     -13.321760917864202,      # [1]: σ_22
            #     -13.321760917864202,      # [2]: σ_33
            #     1.2074287379946301e-8,    # [3]: σ_23
            #     1.2074287379946301e-8,    # [4]: σ_13
            #     1.2074287379946301e-8,    # [5]: σ_12
            # ]
            # [26.643521871951265,    -13.321760917864202,   -13.321760917864202,   1.2074287379946301e-8, 1.2074287379946301e-8, 1.2074287379946301e-8,]
            dj2 = 0.5((sum(ds[[1, 4, 6]] .^ 2.0))
                + 2.0(sum(ds[[2, 3, 5]] .^ 2.0)))
            dj3 =   (ds[1]*(ds[4]*ds[6]-ds[3]*ds[3])
                  - ds[2]*(ds[2]*ds[6]-ds[3]*ds[5])
                  + ds[5]*(ds[2]*ds[3]-ds[4]*ds[5]))
            # di1 = I₁(σ̲̲)
            # dj2 = I₂(ds)
            # dj3 = I₃(ds)
            # @show ds, di1, dj2, dj3
            JJ1 = (dj3^2.0) / (dj2^3.0)
            JJ2 = (dj3    ) / (dj2^1.5)
            JJ3 = (di1    ) / (dj2^0.5)
            # [20250401T1042] (JMA3): this comment (v) is from HEC's original code
            # here I controlled stress triaxiality to 1 (tension (Horstemeyer et al., 2000))
            JJ3 = 1.0
        ##--- nucleation (RK4 integration)
            ddff= ( 𝒹 ^ 0.5 )  /  ( 𝒻 ^ (1.0/3.0) )
            Δη₀ = ϵ̲̲̇′_mag   *   ddff   /   Kic   *   (  a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (
                c * damirr * abs(JJ3) )  )   *   exp(  Tnuc  /  θ  )
                #+ pcc*(1.+sinh(kp1*Si))*abs(JJ3))*exp(pTnuc/θ)
            k₁  = Δη₀  *    η
            k₂  = Δη₀  *  ( η + (0.5k₁*Δt) )
            k₃  = Δη₀  *  ( η + (0.5k₂*Δt) )
            k₄  = Δη₀  *  ( η + (   k₃*Δt) )
            η  += 6.0  \  Δt  *  ( k₁ + 2.0(k₂+k₃) + k₄ ) # ! update ISV
            Δη  = η   *   ϵ̲̲̇′_mag   *   ddff   /   Kic   *   (  a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (
                c * damirr * abs(JJ3) )  )   *   exp(  Tnuc  /  θ  )

            ### Implementation (Horstemeyer et al., 2000)
            #Nuc[i] = pCnuc*exp(TEm[i-1]*ddff/pKic*(paa*(4./27.-JJ1) + pbb*(JJ2) \
            #       + pcc*damirr*abs(JJ3))*exp(pTnuc/θ))

            #nuc0 = PE[i]*ddff/pKic*(paa*(4./27.-JJ1) + pbb*(JJ2) \
            #     + pcc*(1.+sinh(pccsi*Si))*abs(JJ3))*exp(pTnuc/θ)
            #Nuc[i] = pCnuc*exp(nuc0)
        ##--- growth
            ### Implementation (Euler method)
            #dvod = 4./3.*((sqrt(3.)/2.*prr0*ddd/(1.-pnn) \
            #     * sinh(sqrt(3.)*(1.-pnn)*sqrt(2.)/3.*JJ3)) \
            #     * exp(pTgrw*θ))^3              
            #Vod[i] = Vod[i-1] + dvod*dt

            ### Implementation (Horstemeyer et al., 2000)
            νᵥ₀ = νᵥ
            # ! update ISV
            νᵥ = (    4.0    /    3.0    )     *     (#={=#    (   R₀   *   exp(#=[=#
                    ϵ̲̲′_mag  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-nn) )  *  sinh(
                        sqrt(3.0) * (1.0-nn) * sqrt(2.0) / 3.0 * JJ3 )  *  exp( Tgrw * θ )
                #=]=#)   )    ^    3.0    #=}=#)
            Δνᵥ = νᵥ - νᵥ₀

            #vod0 = PE[i]*sqrt(3.)/(2.*(1.-pnn)) \
            #     * sinh(sqrt(3.)*(1.-pnn)*sqrt(2.)/3.*JJ3)*exp(pTgrw*θ)
            #Vod[i] = 4./3.*(prr0*exp(vod0))^3
        ##--- coalesence
        C = 1.0 # ! update ISV
        ##--- damage rate
        ϕ̇ = (Δη*νᵥ) + (η*Δνᵥ) # ! update ISV
        ##--- total damage at current step
        # ? [2025T1048] (JMA3): why the blazes does this phi have 3 re-assignments?
        ϕ₀    = ϕ
        ϕ    += ϕ̇*Δt
        ϕ     = C*η*νᵥ
        ϕ     = max( min(ϕ,0.99999) , 0.0000000001 ) # ! update ISV

        ϕ̇ = (ϕ-ϕ₀) / Δt # ! update ISV

        ϵ̲̲⁽ᴴ⁾ = JJ3

        # ? [20250401T1510] (JMA3): I'm not really sure what this section does either
        ## final process for returning variables
            # TEm[i] = TE[0][i]^2 + TE[1][i]^2 + TE[2][i]^2 \
            #       +(TE[3][i]^2 + TE[4][i]^2 + TE[5][i]^2)*2.
            # TEm[i] = sqrt(TEm[i])*sqrt_twothirds
            ϵ̲̲′_mag += ϵ̲̲̇′_mag*Δt
    
            # Alm[i] = Al[0][i]^2 + Al[1][i]^2 + Al[2][i]^2 \
            #         +(Al[3][i]^2 + Al[4][i]^2 + Al[5][i]^2)*2.
            # Alm[i] = sqrt(Alm[i])*sqrt_threehalves
            α̲̲_mag = sqrt_threehalves * norm_symvec(α̲̲)
        ## saturation stress
            # ϵ̲̲̇′_mag
            ϵ̲̲̲̇′⁽ᵖ⁾_mag = Δγ / Δt * sqrt_twothirds
            #ddp = ddd
            #Alsat = Be + Yp + sqrt(h*ddp/(rd*ddp+rs))
            α̲̲ₛₐₜ_mag = 0.0
            κₛₐₜ  = (  ( (1.0-X) ^ NK )  *  H  *  Hir  *  ϵ̲̲̲̇′⁽ᵖ⁾_mag  /  (
                (sqrt_twothirds*Rdc*ϵ̲̲̲̇′⁽ᵖ⁾_mag) + Rs )  )   ^   (  1.0  /  NK  )
            # if(i >= incnum0-1): Ksat  = vM[i]
            vMₛₐₜ = Be + Y + Yₚ + α̲̲ₛₐₜ_mag + κₛₐₜ
            vMₛₐₜ = κₛₐₜ + Be
    end
    # return (vM,ϵ̲̲′_mag,α̲̲_mag,κ,X,d,ϕ,η,νᵥ,vMₛₐₜ,ϵ̲̲̲̇′⁽ᵖ⁾_mag,t)
    return σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, κₛ, ϕ, η, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d
end

function ContinuumMechanicsBase.predict(
            ψ   ::Cho2019Unified{T}, # , S},
            test::BammannChiesaJohnsonPlasticity.AbstractBCJMetalTest{T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat} # , S<:SymmetricTensor{2, 3, T}}
    M = ψ.N + 1
    # irradiation before damage
    Tirr = ψ.P
    M0, Si, damirr = 0.0, 0.0, 1.0
    if Tirr != 0.0
        kr = kr1 * exp(krt/Tirr)
        Si = (kr*flu) ^ (1.0/kr2)
        M0 = Kr3 * Si
        damirr = exp(  ( kp1 * exp(kpt/Tirr) * flu )  ^  ( 1.0 / kp2 )  )
    end
    # observable state variables
    σ̲̲       = zeros(T, 6)   # Cauchy stress
    ϵ̲̲       = zeros(T, 6)   # total strain
    # internal state variables
    ϵ̲̲⁽ᵖ⁾    = zeros(T, 6)   # plastic strain
    α̲̲       = fill(1e-7, 6) # kinematic hardening
    κ       = M0            # isotropic hardening
    κₛ      = M0            # irradiation hardening
    ## damage
    ϕ       = 1.0e-5        # damage
    η       = ψ.η₀          # void nucleation
    νᵥ      = 0.0           # void growth
    ϕ̇       = 1.0e-5        # damage rate
    ## recrystallization
    X       = 1.0e-10       # total dislocation-free volume fraction
    XR      = 0.0           # total recrystallized volume fraction
    XH      = 0.0           # total reduction of recrystallized volume fraction
    Xd      = 0.5e-10       # total dynamically recrystallized volume fraction
    Xs      = 0.5e-10       # total statically recrystallized volume fraction
    d       = ψ.d₀          # average grain size

    # begin prediction
    σ⃗ = []; push!(σ⃗, σ̲̲)
    ϵ⃗ = []; push!(ϵ⃗, ϵ̲̲)
    t       = 0.0
    for i ∈ range(2, M)
        t += ψ.Δt
        #                                                         update(ψ, t, σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, κₛ, Si, ϕ, η, damirr, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d, (;
        σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, κₛ, ϕ, η, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d = update(ψ, t, σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, κₛ, Si, ϕ, damirr, η, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d, p)
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
    # BCJ-plasticity
    ## yield surface
    # base, exponent
    :C₁,     :C₂,             # V
    :C₃,     :C₄,             # Y
    :C₅,     :C₆,             # f
    ## pressure-dependent yield surface
    :Pₖ₁, :Pₖ₂, :Pₖ₃,
    ## kinematic hardening
    # base, exponent, pressure
    :C₇,     :C₈,     :C₂₁,    # r_d
    :C₉,     :C₁₀,    :C₂₂,    # h
    :C₁₁,    :C₁₂,    :C₂₃,    # r_s
    ## isotropic hardening
    # base, exponent, pressure
    :C₁₃,    :C₁₄,    :C₂₄,    # R_d
    :C₁₅,    :C₁₆,    :C₂₅,    # H
    :C₁₇,    :C₁₈,    :C₂₆,    # R_s
                    :NK,     # * [20250402T1521] (JMA3): I think this is the modifier for finding the k-root
                            # *                         (see Eq. 4.22 in HEC dissertation)
                            # *                         (c. f. `optimize.py` that NK=2.0 by default)
    ## torsion, tension/compression
    :ca, :cb,
    ## dynamic recrystallization
    :Cx1, :Cx2, :Cdp,
    :Cx3, :Cx4, :Csp,
    :Cx5, :Cxa, :Cxb, :Cxc,
    ## static RX (grain growth)
    :n, :ω₀, # E⁺, V⁺, R,
    ## grain size
    # d₀, Cg1, Cg2, Cg3, z,
    :Cg1, :Cg2, :Cg3, :z,
    ## damage
    ### nucleation
    # 𝒹, 𝒻, Kic, a, b, c,
    :a, :b, :c,
    # Cnuc, Tnuc, R₀, nn, Tgrw,
    :pCnuc, :Tnuc, :nn, :Tgrw,
    ## irradiation hardening
    :kr1, :krt, :kr2, :kr3, :kp1, :kpt, :kp2
)

nothing
# end # end of module
