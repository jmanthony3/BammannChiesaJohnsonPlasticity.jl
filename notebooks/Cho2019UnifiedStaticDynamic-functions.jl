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
function update(ψ::Cho2019Unified, σ̲̲, α̲̲, κ, κₛ, ϕ, X, Xd, Xs, XR, XH, d, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, t, (;
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
        # μ       = ψ.μ
        ϵₙ      = ψ.ϵₙ
        N       = ψ.N
        Δϵ      = ψ.Δϵ
        Δt      = ψ.Δt
        M       = N + 1
        T       = typeof(float(θ))
    # calculation constants/functions
        sqrt_twothirds = √(2.0/3.0)
        sqrt_threehalves = √(3.0/2.0)
        p(σ)= sum(σ[[1, 4, 6]]) / 3.0
    # irradiation before damage
        Tirr = pres
        M0, Si, damirr = 0.0, 0.0, 1.0
        if Tirr != 0.0
            kr = kr1 * exp( krt/Tirr )
            Si = (kr*flu) ^ (1.0/kr2)
            M0 = Kr3 * Si
            damirr = exp(  ( kp1 * exp(kpt/Tirr) * flu )  ^  ( 1.0 / kp2 )  )
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
            KT0   = K0      + (            dK0dT * (θ-ttop)    )
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
                # Update Solution
                RRT0 -= FF / dF
                # Convergence Check
                err = abs(dRRT0)
                err <= convg ? break : Niter += 1
                Niter >= (itmax - 1) ? println("BM convergence issue! ", err) : nothing
            end
            # pressure-temperature dependent reference density
            ρ = RRT0 * RT0
    # shear modulus
        #--- 3rd-order Finite Strain (Birch-Murnaghan EOS)
        KT0 = K0     +   (          dK0dT * (θ-ttop)    )
        RT0 = rho0   *   (  1.0 - (   alp * (θ-ttop) )  )
        GT0 = G0     +   (          dG0dT * (θ-ttop)    )
        b1  = (3KT0*dG0dP) - 5GT0
        b2  = 9.0(   (  KT0  ^  2.0  )   *   (  ddG0ddP  +  (
                (1.0/KT0) * (dKdP-4.0) * dG0dP )  )   +    (  35.0GT0  /  9.0  )   )
        f   = 0.5(  ( (ρ/RT0) ^ (2.0/3.0) )  -  1.0  )
        μ   = max(  0.01,  ( (1.0+2.0f) ^ 2.5 )  *  ( GT0 + (b1*f) + 0.5b2 * (f^2.0) )  )
        ν   = 0.3
        G   = (2.0/3.0) * μ * (1.0+ν) / (1.0-2ν)
        if     imat == 1    # OFHC Cu (irradation-ISV model)
           μ = 5.47e4    - (34.1*θ)
           G = 70000.0
        elseif imat == 2    # T91 ferritic steel (Barrett et al., 2018)
           μ = 1.01e5    - (65.0*θ)
           G = 170000.0
        elseif imat == 3    # Ti6Al4V (Hukuhara&Sanpei,1993)
           μ = 4.5e4     - (20.0*θ)
           G = 85000.0
        end
    # deviatoric strain and effective strain rate
        davg = p(Δϵ)
        ddd = sqrt_twothirds * norm_symvec(ϵ̇_eff) / Δt
    # trial damage
        dam1 = 1.0 - ϕ
        dam2 = 1.0 - min(1.0, ϕ̇*Δt/dam1)
    # hydrostatic pressure
        if pres > 0.0
            P = pres
        else
            P_H = (p(σ̲̲)*dam2) + (3.0G*davg*dam1)
        end
    # deviatoric stress and invariants
        di1 = I₁(σ̲̲)
        dj2 = I₂(σ̲̲)
        dj3 = I₃(σ̲̲)
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
        Rdc = Rd * (ddd^-0.0)
    # yield surface parameters
        #iYS: 0-Pressure sensitive (Shear-Mises);
        #     1-Pressure insensitive (Mises);
        #     2-Pressure sensitive (TANH)
        if iYS == 0
            Pa = Pk1   *   (  ( 1.0 + exp(-Pk2/θ) )  ^  ( -Pk3 )  )
            Pc = 0.0Pa
            Pd = Pa - Pc
            imode, Yp = if P <= Pa
                imode, Ft = (P<Pc)   ?   (1, 0.0)   :   (
                    2,  ( 0.5Pd )  *  ( (P-Pc) ^ 2.0 )  *  tanB   )
                Yp = (P*tanB) - Ft
                (   imode,                             Yp   )
            else
                (       3,   (  Pa  -  0.5Pd  )   *  tanB   )
            end
        elseif iYS == 1
            Yp = 0.
        elseif iYS == 2
            #Yp = Pk1*exp(-Pk2*θ)*tanh(B2*P[i])
            #Yp = (Pk1*(1. + exp(-Pk2/θ))^(-Pk3))*tanh(B2*P[i])
            B1 = max(  1.e-10,    ( Pk1 - (Pk2*θ) )  )
            Yp = B1  *  tanh( (Pk3/B1) * P )
        else
            error("iYS > 2 which is not supported.")
        end
    # viscous stress
        #Be = V*log((ddd + sqrt(ddd^2 + F^2))/F)
        #... Using sinh^n for strain rate-stress curve's smooth connection
        Be = V      *      log(     F     \     (#={=#
                (   ddd   ^   (  1.0  /  NK  )   )    +    sqrt(#=[=#
                    (  ( ddd ^ (1.0/NK) )  ^  2.0  )  +  (  F  ^  2.0  )   #=]=#)
            #=}=#)     )
    # previous alpha magnitude
        α_mag   = norm_symvec(α̲̲)
        α_mag  *= sqrt_threehalves
    # REX Model
        ## REX calculation: separated DRX and SRX equations
            if     iREXmethod == 0 # Euler Method (explicit)
                KAlMu   = μ  \  ( (κ^2.0) + (α_mag^2.0) )
                dAlpha  = (  h  *  ddd  )   -   (  ( (sqrt_twothirds* rd*ddd) + rs )  *  ( α_mag ^ 2.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  ddd  )   -   (  ( (sqrt_twothirds*Rdc*ddd) + Rs )  *  (     κ ^ 2.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ  \  ( dKappa + dAlpha )
                Cd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ddd  *  dt  )
                Cs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  dt  )
                Ch      = Cx5 * KAlMu1 * dt
                # pX0     = Cd + Cs
                # dXR     = pX0*(X[i-1]^Cxa)*(1. - X[i-1])^Cxb
                # dXH     = Ch*X[i-1]^Cxc
                # dX      = dXR - dXH
                # new trial
                dXd     = Cd  *  ( X ^ Cxa )  *  ( (1.0-X) ^ Cxb )
                dXs     = Cs  *  ( X ^ Cxa )  *  ( (1.0-X) ^ Cxb )
                dXR     = dXd + dXs
                dXH     = Ch  *  ( X ^ Cxc )
                dX      = dXR - dXH
                Xd     += dXd
                Xs     += dXs
                XH     += dXH
                # dX      = dXR - dXH
                xx      = X + dX
            elseif iREXmethod == 1 # explicit exponential integration algorithm
                KAlMu   = μ  \  ( (κ^2.0) + (α_mag^2.0) )
                dAlpha  = (  h  *  ddd  )   -   (  ( (sqrt_twothirds* rd*ddd) + rs )  *  ( α_mag ^ 2.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  ddd  )   -   (  ( (sqrt_twothirds*Rdc*ddd) + Rs )  *  (     κ ^ 2.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ  \  ( dKappa + dAlpha )
                Cd      = Cx1   *   exp(  -( Cx2 + (   P*Cdp) )  /        θ  )   *   (  KAlMu  *  ddd  *  dt  )
                Cs      = Cx3   *   exp(  -( Cx4 + (   P*Csp) )  /        θ  )   *   (  KAlMu          *  dt  )
                Ch      = Cx5 * KAlMu1 * dt
                Udt     = ( Cd + Cs )  *  ( X ^ Cxa )  *  ( (1.0-X) ^ (Cxb-1.0) )
                Udt     = Udt   +   (  Ch  *  (  X ^ (Cxc-1.0)  )  )
                Vdt     = ( Cd + Cs )  *  ( X ^ Cxa )  *  ( (1.0-X) ^ (Cxb-1.0) )
                xx      = (  X   *   exp( -Udt )  )   +   (  Vdt  *  ( (1.0-exp(-Udt)) / Udt )  )
                pX0     = Cd + Cs
                dXR     = pX0          *  ( X ^ Cxa )  *  ( (1.0-X) ^  Cxb  )
                dXH     = Ch * (X^Cxc)
            elseif iREXmethod == 2 # RK4-explicit method
                # κ = 10.
                KAlMu   = μ  \  ( (κ^1.0) + (α_mag^1.0) )
                dAlpha  = (  h  *  ddd  )   -   (  ( (sqrt_twothirds* rd*ddd) + rs )  *  ( α_mag ^ 3.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  ddd  )   -   (  ( (sqrt_twothirds*Rdc*ddd) + Rs )  *  (     κ ^ 3.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ  \  ( dKappa + dAlpha )
                Cd      = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /  ( R * θ )  )   *   (  KAlMu  *  ddd  )
                # Cd     = Cx1*exp(-(Cx2 + P[i]*Cdp)/θ)*KAlMu*ddd
                Cs      = Cx3   *   exp(  -( Cx4 + (   P*Csp) )  /        θ    )   *      KAlMu
                Ch      = Cx5 * KAlMu1
                CC      = 0.5(Cd + Cs)
                k1      = CC   *   (  X ^ Cxa  )   *   (  ( 1.0-X ) ^ Cxb  )
                k1      = k1 - 0.5(   Ch   *   (    X                 ^  Cxc  )   )
                k2      = CC    *    (   (  X  +     ( dt * k1 )  )   ^   Cxa   )    *    (   (  1.0  -  ( X +    (dt*k1) )  )   ^   Cxb   )
                k2      = k2 - 0.5(   Ch   *   (  ( X +    (dt*k1) )  ^  Cxc  )   )
                k3      = CC    *    (   (  X  +     ( dt * k2 )  )   ^   Cxa   )    *    (   (  1.0  -  ( X +    (dt*k2) )  )   ^   Cxb   )
                k3      = k3 - 0.5(   Ch   *   (  ( X +    (dt*k2) )  ^  Cxc  )   )
                k4      = CC    *    (   (  X  +  2.0( dt * k3 )  )   ^   Cxa   )    *    (   (  1.0  -  ( X + 2.0(dt*k3) )  )   ^   Cxb   )
                k4      = k4 - 0.5(   Ch   *   (  ( X + 2.0(dt*k3) )  ^  Cxc  )   )
                xx      = X   +   (  1.0  /  3.0  )   *   (  ( k1 + 2.0(k2+k3) + k4 )  *  dt  )
                Cd     *= dt
                Cs     *= dt
                Ch     *= dt
                pX0     = Cd + Cs
                dXd     = Cd  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXs     = Cs  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXR     = pX0 *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXH     = Ch * (xx ^ Cxc)
            elseif iREXmethod >= 3 # implicitly solve functions using Newton-Rapson method
                Nitmax = 20
                Ntol   = 1.e-6
                xx     = 0.5
                for k in range(0, Nitmax)
                    if     iREXmethod == 3 # Euler Method (implicit)
                        KAlMu   = μ  \  ( (κ^2.0) + (α_mag^2.0) )
                        dAlpha  = (  h  *  ddd  )   -   (  ( (sqrt_twothirds* rd*ddd) + rs )  *   ( α_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (  H  *  ddd  )   -   (  ( (sqrt_twothirds*Rdc*ddd) + Rs )  *   (     κ ^ 2.0 )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ  \  ( dKappa + dAlpha )
                        Cd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ddd  *  dt  )
                        Cs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  dt  )
                        Ch      = Cx5 * KAlMu1 * dt
                        # TODO [20250331T1119] (JMA3): come back to decrement this section instead
                        f       = X  +   (  ( Cd + Cs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )  )
                        f       = f  - ( Ch * (xx^Cxc) ) - xx
                        df      =    - ( Cd + Cs )  *  Cxb  *  ( (1.0-xx) ^ (Cxb-1.0) )  *  ( xx ^  Cxa )
                        df      = df + ( Cd + Cs )  *  Cxa  *  ( (1.0-xx) ^  Cxb      )  *  ( xx ^ (Cxa-1.0))
                        df      = df - (      Ch    *  Cxc  *  ( xx ^ (Cxc-1.0) )  )
                        df      = df - 1.0
                    elseif iREXmethod == 4 # exponential integration algorithm (asymptotic)
                        KAlMu   = μ  \  ( (κ^2.0) + (α_mag^2.0) )
                        dAlpha  = (  h  *  ddd  )   -   (  ( (sqrt_twothirds* rd*ddd) + rs )  *  ( α_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (  H  *  ddd  )   -   (  ( (sqrt_twothirds*Rdc*ddd) + Rs )  *  (     κ ^ 2.0 )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ  \  ( dKappa + dAlpha )
                        Cd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ddd  *  dt  )
                        Cs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  dt  )
                        Ch      = Cx5 * KAlMu1 * dt
                        Udt     = ( Cd + Cs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        Udt     = Udt    +    (  Ch  *  ( xx ^ (Cxc-1.0) )  )
                        Vdt     = ( Cd + Cs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        f       = (   X * exp(-Udt)   )    +    (   Vdt   *   (  ( 1.0-exp(-Udt) )  /  Udt  )   )    -    xx
                        dUdt    = ( Cxc - 1.0 )  *  Ch  *  ( xx ^ (Cxc-2.0) )
                        dUdt    = dUdt   +   (  ( Cd + Cs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dUdt    = dUdt   -   (  ( Cd + Cs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        dVdt    =            (  ( Cd + Cs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dVdt    = dVdt   -   (  ( Cd + Cs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        df      = -X * dUdt * exp(-Udt)
                        df      = df    +    (   (  ( dVdt/Udt )  -  ( (Vdt*dUdt) / (Udt^2.0) )  )   *   (   1.0  -  exp(-Udt)  )   )
                        df      = df    +    (   (                      Vdt       /  Udt         )   *   (  dUdt  *  exp(-Udt)  )   )
                        df     -= 1.0
                    elseif iREXmethod == 5 # exponential integration algorithm (trapezoidal)
                        KAlMu   = μ  \  ( (κ^2.0) + (α_mag^2.0) )
                        dAlpha  = (  h  *  ddd  )   -   (  ( (sqrt_twothirds* rd*ddd) + rs )  *  ( α_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (  H  *  ddd  )   -   (  ( (sqrt_twothirds*Rdc*ddd) + Rs )  *  (     κ ^ 2.0 )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ  \  ( dKappa + dAlpha )
                        Cd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  θ  )   *   (  KAlMu  *  ddd  *  dt  )
                        Cs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  θ  )   *   (  KAlMu          *  dt  )
                        Ch      = Cx5 * KAlMu1 * dt
                        Udt     = ( Cd + Cs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        Udt     = Udt   +   (  Ch  *  ( xx ^ (Cxc-1.0) )  )
                        U0dt    = ( Cd + Cs )  *  (  X ^ Cxa )  *  ( (1.0- X) ^ (Cxb-1.0) )
                        U0dt    = U0dt  +   (  Ch  *  (  X ^ (Cxc-1.0) )  )
                        Vdt     = ( Cd + Cs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ (Cxb-1.0) )
                        V0dt    = ( Cd + Cs )  *  (  X ^ Cxa )  *  ( (1.0- X) ^ (Cxb-1.0) )
                        f       =                   X  *  exp( -0.5(U0dt+Udt) )
                        f       =  f   +   0.5(  V0dt  *  exp( -0.5(U0dt+Udt) )  +  Vdt  )   -   xx
                        dUdt    = ( Cxc - 1.0 )  *  Ch  *  ( xx ^ (Cxc-2.0) )
                        dUdt    = dUdt   +   (  ( Cd + Cs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dUdt    = dUdt   -   (  ( Cd + Cs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        dVdt    =            (  ( Cd + Cs )  *    Cxa          *  ( xx ^ (Cxa-1.0) )  *  ( (1.0-xx) ^ (Cxb-1.0) )  )
                        dVdt    = dVdt   -   (  ( Cd + Cs )  *  ( Cxb - 1.0 )  *  ( xx ^  Cxa      )  *  ( (1.0-xx) ^ (Cxb-2.0) )  )
                        df      = -0.5dUdt  *  X  *  exp( -0.5(U0dt+Udt) )
                        df      = df   +   0.5(  -0.5dUdt  *  V0dt  *  exp( -0.5(U0dt+Udt) )  +  dVdt  )
                        df      = df   -   1.0
                    end

                    dxx = -f / df
                    xx  = xx + dxx
                    xx  = max(1e-6, min(0.9999999, xx))

                    pX0    = Cd + Cs
                    dXR    = pX0  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                    dXH    =  Ch  *  ( xx ^ Cxc )

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
            X   = xx
            (dX < 0.0)  ?  (X0 = 1.0)  :  nothing # ∵ REX is irreversible
            XR += dXR
            XH += dXH
            #dXd   = dXR*(Cd/pX0)
            #dXs   = dXR*(Cs/pX0)
            Xd += dXd
            Xs += dXs
        ## Grain size kinetics (SGG and grain refinement rate)
            # Grain size rate integration method
            # 0-explicit; 1-implicit; 2-analytic; 3-earlier model (IJP,2019)
            iGSmethod = 0
            di1 = d
            if     iGSmethod == 0 # Forward Euler (explicit)
                # static grain growth rate
                dr      = di1
                dsgk    =  sxk   *   exp(  -( sxE + (1e6P*sxV) )  /  ( R * θ )  )
                dsgg    = dsgk   /   (  sxn  *  ( dr ^ (sxn-1.0) )  )
                # dynamic grain size reduction rate (new version: EPSL2020)
                dred    = Cg1 * X * ddd * (dr^Cg2)
                # total grain size change rate
                d       = dr  +  ( (dsgg-dred) * dt )
                # Z       = ddd*exp((sxE + P[i]*1.e6*sxV)/(R*θ))
                # dss     = (sxk/(Cg3*sxn*0.3))^(1./(sxn-1.+Cg2))*Z^(-(1./(sxn-1.+Cg2)))
            elseif iGSmethod == 1 # Backward Euler: a = 1 (implicit); a = 0.5 (Crank-Nicholson)
                a       = 1.0
                Nitmax  = 20
                Convg   = 1.e-6
                dr      = di1
                # dsgk    = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(R*θ))
                # time downscaling factor for matching to n=4
                tscl    = t  ^  ( (sxn/4.0) - 1.0 )
                dsgk    = sxk   *   exp(  -( sxE + (1e6P*sxV) )  /  ( R * θ )  )   *   tscl
                xx      = dr
                for k in range(0, Nitmax)
                    f   = dr      +      (#={=#     dsgk     *     dt     /     (#=[=#    sxn    *    (
                            (  ( 1.0 - a )  *  ( dr ^ (sxn-1.0) )  )   +   (  a  *  ( xx ^ (sxn-1.0) )  )
                        )    #=]=#)     #=}=#)
                    f  -= Cg1   *   X   *   ddd   *   dt   *   (
                            ( (1.0-a) * (dr^Cg2) )  +  ( a * (xx^Cg2) )  )
                    f  -= xx
                    df  = (  ( dsgk * dt * a * (1.0-sxn) / sxn )  *  ( xx ^ -sxn )  )   -   1.0
                    df -= Cg1    *    X    *    ddd    *    dt    *    (
                            (  Cg2 * a * ( xx ^ (Cg2-1.0) )  )   )
                    dxx = -f / df
                    xx += dxx

                    if abs(dxx) <= Convg
                        break
                    end
                    if k >= Nitmax - 1
                        println("N-R Convg Issue for Grain Size: k >= Nitmax", dxx, k)
                    end
                end
                d       = xx
                prefct  = ( sxk * tscl / (Cg1*sxn*X) )  ^  ( 1.0 / (sxn-1.0+Cg2) )
                dsss    = prefct     *     (    (#=[=#
                        ddd   *   exp(  ( sxE + (1e6P*sxV) ) / ( R * θ )  )
                    #=]=#)    ^    (   -1.0   /   (  sxn  -  1.0  +  Cg2  )   )    )
            elseif iGSmethod == 2 # analytical solution
                # static grain growth
                dsgk    = sxk   *   exp(  -( sxE + (1e6P*sxV) )  /  ( R * θ )  )
                d       = d0    +    (   dsgk   *   t   *   (  t  ^  ( (sxn/4.0) - 1.0 )  )   )    ^    (   1.0   /   sxn   )
            elseif iGSmethod == 3 # original version of DRX grain size kinetics model
                P1      = 300.0
                P2      = 0.18
                P3      = 2.0
                dr      = di1
                tscl    = t  ^  ( (sxn/4.0) - 1.0 )
                dsgk    = sxk   *   exp(  -( sxE + (1e6P*sxV) )  /  ( R * θ )  )   *   tscl
                dssmax  = ( (dsgk*dt) + (dr^sxn) )  ^  ( 1.0 / sxn )
                # ? [20250331T1347] (JMA3): what even is this if-statement?
                if ddd * dt == 0.0
                    dssr = ( (dsgk*dt) + (dr^sxn) )  ^  ( 1.0 / sxn )
                    dssr = dr
                else
                    dss0 = ddd   *   exp(  ( sxE + (1e6P*sxV) )  /  ( R * θ )  )
                    dssr = P1 * (dss0^-P2)
                end
                # ? [20250331T1350] (JMA3): why the addition, subtraction, and increment?
                ddgrw   = (  ( (dsgk*dt) + (dr^sxn) )  ^  (1.0/sxn)  )   -   dr
                dr     += dr + ddgrw
                dss     = min(dssr, dr)
                ddred   = -P3 * X * ddd * dt * dr * (dr-dss)
                # ddred = -Cg3*Xd[i]*dr*(dr - dss)*ddd*dt
                ddred   = max((dss-dr), ddred)
                d       = dr + ddred
            else
                error("iGSmethod > 3 not supported")
            end
        ## Hall-Petch effect
            idzz = 0
            dzz1, dzz0 = if idzz == 0
                ( (d0/d) ^ zz,            1.0 )
            elseif idzz == 1
                (         1.0,   (di1/d) ^ zz )
            elseif idzz == 2
                ( (d0/d)     ,   (di1/d)      ) .^ zz
            else
                error("idzz > 2 which is not supported.")
            end
            # d0 = 1. !Turn on if absolute grain size-stress relation is used
            # YT  = YT*dzz1 ! Turn on if grain size dependent yield is used
    # total time thus far
    t += Δt
    # elastic prediction
        twoμ = 2.0μ
        #--- trial deviatoric stress
            for k in range(0, 3)
                Se[k]  = twoμ*DE[k]
                Str[k] = (S[k][i-1]*dam2) + (Se[k]*dam1)
            end
            for k in range(3, 6)
                Se[k]  = twoμ*DE[k]
                Str[k] = (S[k][i-1]*dam2) + (Se[k]*dam1)
            end
        #--- use of Newton Method for DG and Kappa
        iNewton = 0
        #--- irradiation hardening effect
        Hir = (1.0+Si) ^ 2.0
        #--- trial kappa
        if iNewton == 0
            rdrsk   = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*ddd) )  *  dt  *  ( κ ^ (NK-1.0) )  *  dzz1  )
            Ktr     = κ * X0 * dzz0 / rdrsk
        end
        #--- trial M in isotropic hardening (output only)
        rdrssk  = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*ddd) )  *  dt  *  ( κₛ ^ (NK-1.0) ) *  dzz1  )
        Kstr    = κₛ * X0 * dzz0 / rdrssk
        d       =           (  ( Rs + (sqrt_twothirds*Rdc*ddd) )  *         ( κₛ ^ (NK-1.0) )          )
        if iNewton == 1 # Newton iteration (Backward Euler)
            Nitmax  = 20
            Ntol    = 1.e-06
            Rx0     = X0
            xx      = κ
            for k in range(0, Nitmax)
                RSRD    = 1.0  +  ( (sqrt_twothirds*Rdc*ddd) + Rs ) * dt * ( xx ^ (NK-1.0) )
                F1      = (Rx0*κ/RSRD) - xx
                dF1     = (  -Rx0  *  κ  *  ( (sqrt_twothirds*Rdc*ddd) + Rs )  *  dt  *  ( NK - 1.0 )  *  ( xx ^ (NK-2.0) )  /  ( RSRD ^ 2.0 )  )   -   1.0
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
            Ktr     = xx
            rdrsk   = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*ddd) )  *  dt  *  ( xx ^ (NK-1.0) )  )
        end
        #--- trial alpha
        # * [20250331T1452] (JMA3): this is already defined above
        # Al_mag = Al[0][i-1]^2 + Al[1][i-1]^2 + Al[2][i-1]^2 \
        #        +(Al[3][i-1]^2 + Al[4][i-1]^2 + Al[5][i-1]^2)*2.
        # Al_mag = sqrt(Al_mag)*sqrt_threehalves
        rdrsa   = 1.0   +   (  ( rs + (sqrt_twothirds* rd*ddd) )  *  dt  *  Al_mag            *  dzz1  )
        for k in range(0, 6)
            Altr[k] = Al[k][i-1] * X0 * dzz0 / rdrsa
        end
        #--- Plastic direction tensor N
        for k in range(0, 6)
            Xi[k] = Str[k] - ( (2.0/3.0) * Altr[k] )
        end
        Xi_mag2 = Xi[0]^2 + Xi[1]^2 + Xi[2]^2 \
                +(Xi[3]^2 + Xi[4]^2 + Xi[5]^2)*2.0
        Xi_mag = sqrt(Xi_mag2)
        for k in range(0, 6)
            N[k] = Xi[k] / Xi_mag
        end
    # check plasticity
        ak     = YT + Ktr + Be + Yp
        critra = Xi_mag - (sqrt_twothirds*ak*dam1)
        #if(TEm[i-1] > 0.02): critra = 100.
    # Radial-Return
    if critra <= 0.0 # elastic solution update
        EPflag = 1
        # deviatoric stress update
        for k in range(0, 6)
            S[k][i]   = Str[k]
            Sig[k][i] = S[k][i]
        end
        # Cauchy stress update
        for k in range(0, 3)
            Sig[k][i] = Sig[k][i] + P_H
        end
        # von Mises stress update
        vM[i]  = S[0][i]^2 + S[1][i]^2 + S[2][i]^2 \
                +(S[3][i]^2 + S[4][i]^2 + S[5][i]^2)*2.
        vM[i]  = sqrt(vM[i])*sqrt_threehalves
        # kinematic hardening & total strain update
        for k in range(0, 6)
            Al[k][i] = Altr[k]
            TE[k][i] = TE[k][i-1] + DE[k]
        end
        # kinematic hardening update
        Alm[i] = Al[0][i]^2 + Al[1][i]^2 + Al[2][i]^2 \
                +(Al[3][i]^2 + Al[4][i]^2 + Al[5][i]^2)*2.
        Alm[i] = sqrt(Alm[i])*sqrt_threehalves
        # isotropic hardening update
        K[i]   = Ktr
        # irradiation hardening update
        Ks[i]  = Kstr
        # total equivalent strain update
        TEm[i] = TEm[i-1] + (ddd*dt)
        # volumetric strain update
        VE[i]  = VE[i-1] + (3.0*davg)

        # damage update
        Phi[i] = Phi[i-1]
        Nuc[i] = Nuc[i-1]
        Vod[i] = Vod[i-1]
        dPhi[i] = 0.0

        Alsat = 0.0
        Ksat  = (  (H * Hir * ddd )  /  ( (sqrt_twothirds*Rdc*ddd) + Rs )  )   ^   (  1.0  /  NK  )
        ddp   = ddd
        vMsat = Be + YT + Yp + Alsat + Ksat
    else # plastic solution (Radial return starts)
        EPflag = 2
        #--- Plastic strain increment solution
        if     iNewton == 0 # analytical solution for DG
            DG = (    Xi_mag    -    (   sqrt_twothirds   *   ak   *   dam1   )    )     /     (
                (   dam1   *   twoμ   )    +    (   dam1   *   (  2.0  /  3.0  )   *   (
                        ( (1.0-X[i]) ^ NK )  *  dzz1  *  ( (h/rdrsa) + (H*Hir/rdrsk) )  )   )    )
        elseif iNewton == 1 # Newton-Rapson for DG and Kappa
            DG = (    Xi_mag    -    (   sqrt_twothirds   *   ak   *   dam1   )    )     /     (
                (   dam1   *   twoμ   )    +    (   dam1   *   (  2.0  /  3.0  )   *   (
                        ( (1.0-X[i]) ^ NK )  *  dzz1  *  ( (h/rdrsa) + (H*Hir/rdrsk) )  )   )    )

            Nitmax  = 20
            Ntol    = 1e-6
            xx1     = DG
            xx2     = κ
            K0      = κ
            Rx      = (1.0-X[i]) ^ NK
            th      = 1.0 # 1-Backward Euler; 0.5-Midpoint; 0-Forward Euler
            for k in range(0, Nitmax)
                thK0thK     = (  ( 1.0 - th )  *  ( K0 ^ (NK-1.0) )  )   +   (  th  *  ( xx2 ^ (NK-1.0) )  )
                Rdxx1Rsdt   = 1.0  +  (  ( (sqrt_twothirds*Rdc*xx1) + (Rs*dt) )  *  thK0thK  )
                F1      = Xi_mag   -   (  twoμ * xx1  )   -   (
                        sqrt_twothirds  *  ( K0 + (Rx*H*Hir*xx1) )  /  Rdxx1Rsdt  )   -   (
                        sqrt_twothirds  *  ( Be + YT + Yp )  )
                dF1dx1  = -twoμ     -     (    sqrt_twothirds    *    (
                    (  Rx  *  H  *  Hir  /  Rdxx1Rsdt  )   -   (
                        ( K0 + (Rx*H*Hir*xx1) )  *  sqrt_twothirds  *  Rdc  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )   )    )
                dF1dx2  = -1.0 * sqrt_twothirds
                F2      = (  ( K0 + (Rx*H*Hir*xx1) )  /  Rdxx1Rsdt  )   -   xx2
                dF2dx1  = (  Rx  *  H  *  Hir  /  Rdxx1Rsdt  )   -   (
                    ( K0 + (Rx*H*Hir*xx1) )  *  sqrt_twothirds  *  Rdc  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )
                dF2dx2  = -( K0 + (Rx*H*Hir*xx1) )  *  (
                    (sqrt_twothirds*Rdc*xx1) + (Rs*dt) )  *  th  *  ( NK - 1.0 )  *  ( xx2 ^ (NK-2.0) )
                dF2dx2  = ( dF2dx2 / (Rdxx1Rsdt^2.0) )  -  1.0

                a11    = dF1dx1
                a12    = dF1dx2
                a21    = dF2dx1
                a22    = dF2dx2
                dxx2   = ( (-F2*a11/a21) + F1 )  /  ( (a22 * (a11/a21)) - a12 )
                dxx1   = ( (-a12*dxx2)   - F1 )  /  a11

                xx1n   = xx1
                xx2n   = xx2
                
                xx1    = xx1 + dxx1
                xx2    = xx2 + dxx2

                if abs(dxx1/xx1n) <= Ntol && abs(dxx2/xx2n) <= Ntol
                    break
                end
                if k >= Nitmax - 1
                    println("Gamma-K: N-R Conv. Issue: k >= Nitmax")
                end
            end
            DG   = xx1
            K[i] = xx2
        end
        #--- stress solution
            # deviatoric stress update
            for k in range(0, 6)
                S[k][i] = Str[k] - (dam1 * twoμ * DG * N[k])
            end
            # Cauchy stress update
            for k in range(0, 3)
                Sig[k][i] = S[k][i] + P_H
            end
            for k in range(3, 6)
                Sig[k][i] = S[k][i]
            end
            # von Mises stress update
            vM[i] = S[0][i]^2 + S[1][i]^2 + S[2][i]^2 \
                    +(S[3][i]^2 + S[4][i]^2 + S[5][i]^2)*2.
            vM[i] = sqrt(vM[i])*sqrt_threehalves
        #--- total deviatoric strain
        for k in range(0, 6)
            TE[k][i] = TE[k][i-1] + DE[k]
        end
        #--- total plastic strain
        PE[i] = PE[i-1] + (sqrt_twothirds * DG)
        #--- total volumetric strain
        VE[i] = VE[i-1] + (3.0 * davg)
        #--- alpha solution
        for k in range(0, 6)
            Al[k][i] = Altr[k]   +   (
                ( (1.0-X[i]) ^ NK )  *  dzz1  *  dam1  *  h  *  DG  *  N[k]  /  rdrsa  )
        end
        #--- kappa solution
        K[i] = (   iNewton   !=   0   )    ?    (   xx2   )    :    (#=[=#   Ktr   +   (
            ( (1.0-X[i]) ^ NK )  *  sqrt_twothirds  *  dzz1  *  dam1  *  H  *  Hir  *  DG  /  rdrsk  )   #=]=#)
        #--- irradiation hardening solution in isotropic hardening
        Sir   = Si^2
        Ks[i] = Kstr   +   (
            ( (1.0-X[i]) ^ NK )  *  sqrt_twothirds  *  dzz1  *  dam1  *  H  *  Sir  *  DG  /  rdrssk  )
        #--- damage
            for k in range(0, 3)
                ds[k] = S[k][i]
            end
            for k in range(3, 6)
                ds[k] = S[k][i]
            end
            di1  = Sig[0][i] + Sig[1][i] + Sig[2][i]
            #di1  = 3.0 * P_H
            dj2  = 0.5*(ds[0]^2 + ds[1]^2 + ds[2]^2 \
                + (ds[3]^2 + ds[4]^2 + ds[5]^2)*2.)
            dj3  = ds[0]*(ds[1]*ds[2]-ds[4]*ds[4])       \
                - ds[3]*(ds[3]*ds[2]-ds[4]*ds[5])\
                + ds[5]*(ds[3]*ds[4]-ds[1]*ds[5])
            JJ1  = (dj3^2.0) / (dj2^3.0)
            JJ2  = (dj3)       / (dj2^1.5)
            JJ3  = (di1)       / (dj2^0.5)
            # here I controlled stress triaxiality to 1 (tension (Horstemeyer et al., 2000))
            JJ3  = 1.
        ##--- nucleation (RK4 integration)
            ddff = ( pdd ^ 0.5 )  /  ( pff ^ (1.0/3.0) )
            dnuc0= ddd   *   ddff   /   pKic   *   (  paa  *  ( (4.0/27.0) - JJ1 )  +  ( pbb * JJ2 )  +  (
                pcc * damirr * abs(JJ3) )  )   *   exp(  pTnuc  /  θ  )
                #+ pcc*(1.+sinh(kp1*Si))*abs(JJ3))*exp(pTnuc/θ)
            k1   = dnuc0  *    Nuc[i-1]
            k2   = dnuc0  *  ( Nuc[i-1] + (0.5k1*dt) )
            k3   = dnuc0  *  ( Nuc[i-1] + (0.5k2*dt) )
            k4   = dnuc0  *  ( Nuc[i-1] + (   k3*dt) )
            Nuc[i] = Nuc[i-1]   +   (  1.0  /  6.0  *  dt  *  ( k1 + 2.0(k2+k3) + k4 )  )
            dnuc = Nuc[i]   *   ddd   *   ddff   /   pKic   *   (  paa  *  ( (4.0/27.0) - JJ1 )  +  ( pbb * JJ2 )  +  (
                pcc * damirr * abs(JJ3) )  )   *   exp(  pTnuc  /  θ  )

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
            Vod[i] = (    4.0    /    3.0    )     *     (#={=#    (   prr0   *   exp(#=[=#
                    TEm[i-1]  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-pnn) )  *  sinh(
                        sqrt(3.0) * (1.0-pnn) * sqrt(2.0) / 3.0 * JJ3 )  *  exp( pTgrw * θ )
                #=]=#)   )    ^    3.0    #=}=#)
            dvod = Vod[i] - Vod[i-1]

            #vod0 = PE[i]*sqrt(3.)/(2.*(1.-pnn)) \
            #     * sinh(sqrt(3.)*(1.-pnn)*sqrt(2.)/3.*JJ3)*exp(pTgrw*θ)
            #Vod[i] = 4./3.*(prr0*exp(vod0))^3
        ##--- coalesence
        coal = 1.
        ##--- damage rate
        dPhi[i] = (dnuc*Vod[i]) + (Nuc[i]*dvod)
        ##--- total damage at current step
        Phi[i] = Nuc[i]*Vod[i]*coal
        Phi[i] = Phi[i-1] + (dPhi[i]*dt)
        Phi[i] = max(min(Phi[i],0.99999),0.0000000001)

        dPhi[i] = (Phi[i]-Phi[i-1]) / dt

        VE[i] = JJ3
    end

    ## final process for returning variables
        #TEm[i] = TE[0][i]^2 + TE[1][i]^2 + TE[2][i]^2 \
        #       +(TE[3][i]^2 + TE[4][i]^2 + TE[5][i]^2)*2.
        #TEm[i] = sqrt(TEm[i])*sqrt_twothirds
        TEm[i] = TEm[i-1] + (ddd*dt)

        Alm[i] = Al[0][i]^2 + Al[1][i]^2 + Al[2][i]^2 \
                +(Al[3][i]^2 + Al[4][i]^2 + Al[5][i]^2)*2.
        Alm[i] = sqrt(Alm[i])*sqrt_threehalves

    ## saturation stress
        ddp = DG / dt * sqrt_twothirds
        #ddp = ddd
        #Alsat = Be + Yp + sqrt(h*ddp/(rd*ddp+rs))
        Alsat = 0.0
        Ksat  = (  ( (1.0-X[i]) ^ NK )  *  H  *  Hir  *  ddp  /  (
            (sqrt_twothirds*Rdc*ddp) + Rs )  )   ^   (  1.0  /  NK  )
        #if(i >= incnum0-1): Ksat  = vM[i]
        vMsat = Be + YT + Yp + Alsat + Ksat
        vMsat = Ksat + Be
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
