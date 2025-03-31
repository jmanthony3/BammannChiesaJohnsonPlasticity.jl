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
function update(ψ::Cho2019Unified, σ̲̲, α̲̲, κ, ϕ, X, Xd, Xs, XR, XH, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, (;
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
    # yield surface parameters
        #iYS: 0-Pressure sensitive (Shear-Mises);
        #     1-Pressure insensitive (Mises);
        #     2-Pressure sensitive (TANH)
        if iYS == 0
            Pa = Pk1 * (   (  1.0 + exp(-Pk2 / temp)  ) ^ (  -Pk3  )   )
            Pc = 0.0Pa
            Pd = Pa - Pc
            imode, Yp = if P <= Pa
               imode, Ft = P < Pc ? (1, 0.0) : (2, (1.0 / (2.0Pd)) * ((P - Pc) ^ 2) * tanB)
               Yp = (P * tanB) - Ft
               (imode, Yp)
            else
               (3,    ( Pa - (Pd / 2.0) ) * tanB   )
            end
        elseif iYS == 1
            Yp = 0.
        elseif iYS == 2
            #Yp = Pk1*exp(-Pk2*temp)*tanh(B2*P[i])
            #Yp = (Pk1*(1. + exp(-Pk2/temp))^(-Pk3))*tanh(B2*P[i])
            B1 = max(   1.e-10,    ( Pk1 - (Pk2 * temp) )   )
            Yp = B1 * tanh(   ( Pk3 / B1 ) * P   )
        else
            error("iYS > 2 which is not supported.")
        end
    # viscous stress
        #Be = V*log((ddd + sqrt(ddd^2 + F^2))/F)
        #... Using sinh^n for strain rate-stress curve's smooth connection
        Be = V * log(     F   \   #=[=#(   (  ddd ^ (1.0 / NK)  )  +  sqrt(  ( (ddd ^ (1.0 / NK)) ^ 2 )  +  ( F ^ 2 )  )   )#=]=#     )
    # previous alpha magnitude
    α_mag       = norm_symvec(α̲̲)
    α_mag      *= sqrt_threehalves


    # REX Model
        ## REX calculation: separated DRX and SRX equations
            if     iREXmethod == 0 # Euler Method (explicit)
                KAlMu   = μ \ ( (κ ^ 2.0) + (α_mag ^ 2.0) )
                dAlpha  = (   h * ddd   )    -    (   (  ( sqrt_twothirds * rd  * ddd ) + rs  )   *   (  α_mag ^ 2.0  )   )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (   H * ddd   )    -    (   (  ( sqrt_twothirds * Rdc * ddd ) + Rs  )   *   (  κ ^ 2.0  )   )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ \ (dKappa + dAlpha)
                Cd      = Cx1    *    exp(   -(  Cx2 + (P * Cdp)  )  /  θ   )    *    (   KAlMu * ddd * dt   )
                Cs      = Cx3    *    exp(   -(  Cx4 + (P * Csp)  )  /  θ   )    *    (   KAlMu * dt   )
                Ch      = Cx5 * KAlMu1 * dt
                # pX0     = Cd + Cs
                # dXR     = pX0*(X[i-1]^Cxa)*(1. - X[i-1])^Cxb
                # dXH     = Ch*X[i-1]^Cxc
                # dX      = dXR - dXH
                # new trial
                dXd     = Cd  *  ( X ^ Cxa )  *  ( (1.0 - X) ^ Cxb )
                dXs     = Cs  *  ( X ^ Cxa )  *  ( (1.0 - X) ^ Cxb )
                dXR     = dXd + dXs
                dXH     = Ch  *  ( X ^ Cxc )
                dX      = dXR - dXH
                Xd     += dXd
                Xs     += dXs
                XH     += dXH
                # dX      = dXR - dXH
                xx      = X + dX
            elseif iREXmethod == 1 # explicit exponential integration algorithm
                KAlMu   = μ \ ( (κ ^ 2.0) + (α_mag ^ 2.0) )
                dAlpha  = (   h * ddd   )    -    (   (  ( sqrt_twothirds * rd  * ddd ) + rs  )   *   (  α_mag ^ 2.0  )   )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (   H * ddd   )    -    (   (  ( sqrt_twothirds * Rdc * ddd ) + Rs  )   *   (  κ ^ 2.0  )   )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ \ (dKappa + dAlpha)
                Cd      = Cx1    *    exp(   -(  Cx2 + (P * Cdp)  )  /  θ   )    *    (   KAlMu * ddd * dt   )
                Cs      = Cx3    *    exp(   -(  Cx4 + (P * Csp)  )  /  θ   )    *    (   KAlMu * dt   )
                Ch      = Cx5 * KAlMu1 * dt
                Udt     = (  Cd + Cs  )   *   (  X ^ Cxa  )   *   (  ( 1.0 - X )  ^  ( Cxb - 1.0 )  )
                Udt     = Udt    +    (   Ch  *  (  X ^ (Cxc - 1.0)  )   )
                Vdt     = (  Cd + Cs  )   *   (  X ^ Cxa  )   *   (  ( 1.0 - X )  ^  ( Cxb - 1.0 )  )
                xx      = (   X * exp(-Udt)   )    +    (   Vdt  *  (  ( 1.0 - exp(-Udt) ) / Udt  )   )
                pX0     = Cd + Cs
                dXR     = pX0   *   (  X ^ Cxa  )   *   (  (1.0 - X) ^ Cxb  )
                dXH     = Ch * (X ^ Cxc)
            elseif iREXmethod == 2 # RK4-explicit method
                # κ = 10.
                KAlMu   = μ \ ( (κ ^ 1.0) + (α_mag ^ 1.0) )
                dAlpha  = (   h * ddd   )    -    (   (  ( sqrt_twothirds * rd  * ddd ) + rs  )   *   (  α_mag ^ 3.0  )   )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (   H * ddd   )    -    (   (  ( sqrt_twothirds * Rdc * ddd ) + Rs  )   *   (  κ ^ 3.0  )   )
                dKappa  = max(0.0, dKappa)
                KAlMu1  = μ \ (dKappa + dAlpha)
                Cd      = Cx1    *    exp(   -(  Cx2 + (1e6P * Cdp)  )   /   (  R * θ  )   )    *    (   KAlMu * ddd   )
                # Cd     = Cx1*exp(-(Cx2 + P[i]*Cdp)/temp)*KAlMu*ddd
                Cs      = Cx3    *    exp(   -(  Cx4 + (   P * Csp)  )   /          θ      )    *    (   KAlMu   )
                Ch      = Cx5 * KAlMu1
                CC      = 0.5(Cd + Cs)
                k1      = CC   *   (  X ^ Cxa  )   *   (  ( 1.0 - X ) ^ Cxb  )
                k1      = k1 - 0.5(  Ch * ( X ^ Cxc )  )
                k2      = CC   *   (  ( X +    (dt * k1) ) ^ Cxa  )   *   (  ( 1.0 - (X +    (dt * k1)) ) ^ Cxb  )
                k2      = k2 - 0.5(  Ch * ( X +    (dt * k1) ) ^ Cxc  )
                k3      = CC   *   (  ( X +    (dt * k2) ) ^ Cxa  )   *   (  ( 1.0 - (X +    (dt * k2)) ) ^ Cxb  )
                k3      = k3 - 0.5(  Ch * ( X +    (dt * k2) ) ^ Cxc  )
                k4      = CC   *   (  ( X + 2.0(dt * k3) ) ^ Cxa  )   *   (  ( 1.0 - (X + 2.0(dt * k3)) ) ^ Cxb  )
                k4      = k4 - 0.5(  Ch * ( X + 2.0(dt * k3) ) ^ Cxc  )
                xx      = X   +   (  1.0 / 3.0  )   *   (  ( k1  +  2.0(k2 + k3)  +  k4 )  *  dt  )
                Cd     *= dt
                Cs     *= dt
                Ch     *= dt
                pX0     = Cd + Cs
                dXd     = Cd  *  ( xx ^ Cxa )  *  ( (1.0 - xx) ^ Cxb )
                dXs     = Cs  *  ( xx ^ Cxa )  *  ( (1.0 - xx) ^ Cxb )
                dXR     = pX0 *  ( xx ^ Cxa )  *  ( (1.0 - xx) ^ Cxb )
                dXH     = Ch * (xx ^ Cxc)
            elseif iREXmethod >= 3 # implicitly solve functions using Newton-Rapson method
                Nitmax = 20
                Ntol   = 1.e-6
                xx     = 0.5
                for k in range(0, Nitmax)
                    if     iREXmethod == 3 # Euler Method (implicit)
                        KAlMu   = μ \ ( (κ ^ 2.0) + (α_mag ^ 2.0) )
                        dAlpha  = (   h * ddd   )    -    (   (  ( sqrt_twothirds * rd  * ddd ) + rs  )   *   (  α_mag ^ 2.0  )   )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (   H * ddd   )    -    (   (  ( sqrt_twothirds * Rdc * ddd ) + Rs  )   *   (  κ ^ 2.0  )   )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ \ (dKappa + dAlpha)
                        Cd      = Cx1    *    exp(   -(  Cx2 + (P * Cdp)  )  /  θ   )    *    (   KAlMu * ddd * dt   )
                        Cs      = Cx3    *    exp(   -(  Cx4 + (P * Csp)  )  /  θ   )    *    (   KAlMu * dt   )
                        Ch      = Cx5 * KAlMu1 * dt
                        # TODO [20250331T1119] (JMA3): come back to decrement this section instead
                        f       = X  +   (  ( Cd + Cs )  *  ( xx ^ Cxa )  *  ( (1.0 - xx) ^ Cxb )  )
                        f       = f  - ( Ch * (xx ^ Cxc) ) - xx
                        df      =    - ( Cd + Cs )  *  Cxb  *  ( (1.0 - xx) ^ (Cxb - 1.0) )  *  ( xx ^  Cxa )
                        df      = df + ( Cd + Cs )  *  Cxa  *  ( (1.0 - xx) ^  Cxb        )  *  ( xx ^ (Cxa - 1.0))
                        df      = df - (  Ch   *   Cxc   *   ( xx  ^  (Cxc - 1.0) )  )
                        df      = df - 1.0
                    elseif iREXmethod == 4 # exponential integration algorithm (asymptotic)
                        KAlMu   = μ \ ( (κ ^ 2.0) + (α_mag ^ 2.0) )
                        dAlpha  = (   h * ddd   )    -    (   (  ( sqrt_twothirds * rd  * ddd ) + rs  )   *   (  α_mag ^ 2.0  )   )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (   H * ddd   )    -    (   (  ( sqrt_twothirds * Rdc * ddd ) + Rs  )   *   (  κ ^ 2.0  )   )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ \ (dKappa + dAlpha)
                        Cd      = Cx1    *    exp(   -(  Cx2 + (P * Cdp)  )  /  θ   )    *    (   KAlMu * ddd * dt   )
                        Cs      = Cx3    *    exp(   -(  Cx4 + (P * Csp)  )  /  θ   )    *    (   KAlMu * dt   )
                        Ch      = Cx5 * KAlMu1 * dt
                        Udt     = (  Cd + Cs  )   *   (  xx ^ Cxa  )   *   (  ( 1.0 - xx )  ^  ( Cxb - 1.0 )  )
                        Udt     = Udt    +    (   Ch  *  (  xx ^ (Cxc - 1.0)  )   )
                        Vdt     = (  Cd + Cs  )   *   (  xx ^ Cxa  )   *   (  ( 1.0 - xx )  ^  ( Cxb - 1.0 )  )
                        f       = (   X * exp(-Udt)   )    +    (   Vdt  *  (  ( 1.0 - exp(-Udt) ) / Udt  )   )    -    xx
                        dUdt    = ( Cxc - 1.0 )  *  Ch  *  ( xx ^ (Cxc - 2.0) )
                        dUdt    = dUdt   +   (  (Cd + Cs)   *    Cxa          *   (xx ^ (Cxa - 1.0))   *   ( (1.0 - xx) ^ (Cxb - 1.0) )  )
                        dUdt    = dUdt   -   (  (Cd + Cs)   *   (Cxb - 1.0)   *   (xx ^  Cxa       )   *   ( (1.0 - xx) ^ (Cxb - 2.0) )  )
                        dVdt    =            (  (Cd + Cs)   *    Cxa          *   (xx ^ (Cxa - 1.0))   *   ( (1.0 - xx) ^ (Cxb - 1.0) )  )
                        dVdt    = dVdt   -   (  (Cd + Cs)   *   (Cxb - 1.0)   *   (xx ^  Cxa       )   *   ( (1.0 - xx) ^ (Cxb - 2.0) )  )
                        df      = -X * dUdt * exp(-Udt)
                        df      = df   +   (   (  ( dVdt / Udt )  -  ( (Vdt * dUdt) / (Udt ^ 2.0) )  )   *   (  1.0 - exp(-Udt)  )   )
                        df      = df   +   (                           (Vdt /  Udt)                      *   ( dUdt * exp(-Udt)  )   ) - 1.0
                    elseif iREXmethod == 5 # exponential integration algorithm (trapezoidal)
                        KAlMu   = μ \ ( (κ ^ 2.0) + (α_mag ^ 2.0) )
                        dAlpha  = (   h * ddd   )    -    (   (  ( sqrt_twothirds * rd  * ddd ) + rs  )   *   (  α_mag ^ 2.0  )   )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (   H * ddd   )    -    (   (  ( sqrt_twothirds * Rdc * ddd ) + Rs  )   *   (  κ ^ 2.0  )   )
                        dKappa  = max(0.0, dKappa)
                        KAlMu1  = μ \ (dKappa + dAlpha)
                        Cd      = Cx1    *    exp(   -(  Cx2 + (P * Cdp)  )  /  θ   )    *    (   KAlMu * ddd * dt   )
                        Cs      = Cx3    *    exp(   -(  Cx4 + (P * Csp)  )  /  θ   )    *    (   KAlMu * dt   )
                        Ch      = Cx5 * KAlMu1 * dt
                        Udt     = (  Cd + Cs  )   *   (  xx ^ Cxa  )   *   (  ( 1.0 - xx )  ^  ( Cxb - 1.0 )  )
                        Udt     = Udt    +    (   Ch  *  (  xx ^ (Cxc - 1.0)  )   )
                        U0dt    = (  Cd + Cs  )   *   (   X ^ Cxa  )   *   (  ( 1.0 -  X )  ^  ( Cxb - 1.0 )  )
                        U0dt    = U0dt   +    (   Ch  *  (   X ^ (Cxc - 1.0)  )   )
                        Vdt     = (  Cd + Cs  )   *   (  xx ^ Cxa  )   *   (  ( 1.0 - xx )  ^  ( Cxb - 1.0 )  )
                        V0dt    = (  Cd + Cs  )   *   (   X ^ Cxa  )   *   (  ( 1.0 -  X )  ^  ( Cxb - 1.0 )  )
                        f       =                X * exp( -0.5(U0dt + Udt) )
                        f       =  f   +   0.5(             V0dt * exp( -0.5(U0dt + Udt) )  +   Vdt  )   -   xx
                        dUdt    = ( Cxc - 1.0 )  *  Ch  *  ( xx ^ (Cxc - 2.0) )
                        dUdt    = dUdt   +   (  (Cd + Cs)   *    Cxa          *   (xx ^ (Cxa - 1.0))   *   ( (1.0 - xx) ^ (Cxb - 1.0) )  )
                        dUdt    = dUdt   -   (  (Cd + Cs)   *   (Cxb - 1.0)   *   (xx ^  Cxa       )   *   ( (1.0 - xx) ^ (Cxb - 2.0) )  )
                        dVdt    =            (  (Cd + Cs)   *    Cxa          *   (xx ^ (Cxa - 1.0))   *   ( (1.0 - xx) ^ (Cxb - 1.0) )  )
                        dVdt    = dVdt   -   (  (Cd + Cs)   *   (Cxb - 1.0)   *   (xx ^  Cxa       )   *   ( (1.0 - xx) ^ (Cxb - 2.0) )  )
                        df      = -0.5dUdt * X * exp( -0.5(U0dt + Udt) )
                        df      = df   +   0.5(  -0.5dUdt * V0dt * exp( -0.5(U0dt + Udt) )  +  dVdt  )
                        df      = df   -   1.0
                    else
                        error("iREXMethod > 5 which is not supported.")
                    end

                    dxx = -f / df
                    xx  = xx + dxx
                    xx  = max(1e-6, min(0.9999999, xx))

                    pX0    = Cd + Cs
                    dXR    = pX0  *  ( xx ^ Cxa )  *  ( (1.0 - xx) ^ Cxb )
                    dXH    = Ch   *  ( xx ^ Cxc )

                    if abs(dxx) <= Ntol
                        break
                    end
                    if k >= Nitmax-1
                        println("Newton-Rapson Convergence Issue: k >= Nitmax")
                    end
                end
            end

            # Final solution and estimate volume of DRX and SRX
            xx  = max(1e-6, min(0.9999999, xx))
            dX  = xx - X
            X0  = 1.0  -  ( dX / (1.0 - X) )
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

            # Analytical solution
            if(iGSmethod == 2):
                # Static grain growth
                dsgk  = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(R*temp))
                d[i]  = d[0] + (dsgk*Time[i]*Time[i]**(sxn/4.-1))**(1./sxn)

            # Explicit (forward Euler)
            if(iGSmethod == 0):
                # Static grain growth rate
                dr    = d[i-1]
                dsgk  = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(R*temp))
                dsgg  = dsgk/(sxn*dr**(sxn-1.))
                # Dynamic grain size reduction rate (new version: EPSL2020)
                dred  = Cg1*X[i]*ddd*dr**Cg2
                # Total grain size change rate
                d[i]  = dr + (dsgg - dred)*dt
                #Z     = ddd*exp((sxE + P[i]*1.e6*sxV)/(R*temp))
                #dss   = (sxk/(Cg3*sxn*0.3))**(1./(sxn-1.+Cg2))*Z**(-(1./(sxn-1.+Cg2)))

            # Implicit (backward Euler a=1; Crank-Nicholson a=0.5)
            if(iGSmethod == 1):
                a      = 1.0
                Nitmax = 20
                Convg  = 1.e-6
                dr     = d[i-1]
                #dsgk   = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(R*temp))
                # time downscaling factor for matching to n=4
                tscl   = Time[i]**(sxn/4.-1)
                dsgk   = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(R*temp))*tscl
                xx     = dr
                for k in range(0,Nitmax):
                    f   = dr + dsgk*dt/(sxn*((1.-a)*dr**(sxn-1.) + a*xx**(sxn-1.)))
                    f   = f  - Cg1*X[i]*ddd*dt*((1.-a)*dr**Cg2 + a*xx**Cg2)
                    f   = f  - xx
                    df  = (dsgk*dt*a*(1.-sxn)/sxn)*xx**(-sxn) - 1.
                    df  = df - Cg1*X[i]*ddd*dt*(Cg2*a*xx**(Cg2-1.))
                    dxx = -f/df
                    xx  = xx + dxx

                    if(abs(dxx) <= Convg):
                        break
                    if(k >= Nitmax-1):
                        print("N-R Convg Issue for Grain Size: k >= Nitmax", dxx, k)

                d[i]   = xx
                prefct = (sxk*tscl/(Cg1*sxn*X[i]))**(1./(sxn-1.+Cg2))
                dsss   = prefct*(ddd*exp((sxE+P[i]*1.e6*sxV)/(R*temp)))\
                    **(-1./(sxn-1.+Cg2))

            # Original version of DRX grain size kinetics model
            if(iGSmethod == 3):
                P1      = 300.
                P2      = 0.18
                P3      = 2.
                dr      = d[i-1]
                tscl    = Time[i]**(sxn/4.-1)
                dsgk    = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(R*temp))*tscl
                dssmax  = (dsgk*dt + dr**sxn)**(1./sxn)
                if(ddd*dt == 0.):
                    dssr = (dsgk*dt + dr**sxn)**(1./sxn)
                    dssr = dr
                else:
                    dss0 = ddd*exp((sxE + P[i]*1.e6*sxV)/(R*temp))
                    dssr = P1*(dss0**(-P2))

                ddgrw = (dsgk*dt + dr**sxn)**(1./sxn) - dr
                dr    = dr + ddgrw
                dss   = min(dssr, dr)
                ddred = -P3*X[i]*ddd*dt*dr*(dr - dss)
                #ddred = -Cg3*Xd[i]*dr*(dr - dss)*ddd*dt
                ddred = max((dss - dr), ddred)
                d[i]  = dr + ddred

        ## Hall-Petch effect
            idzz = 0
            if (idzz == 0):
                dzz1 = (d0/d[i])**zz
                dzz0 = 1.
            elif (idzz == 1):
                dzz1 = 1.
                dzz0 = (d[i-1]/d[i])**zz
            elif (idzz == 2):
                dzz1 = (d0/d[i])**zz
                dzz0 = (d[i-1]/d[i])**zz
            #d0 = 1. !Turn on if absolute grain size-stress relation is used
            #YT  = YT*dzz1 ! Turn on if grain size dependent yield is used
    # # trial guesses: ISVs (from recovery) and stress
    # recovery    = Δt * (rd * ( sqrt_twothirds * ϵ̇_eff ) + rs) * α_mag     # recovery for alpha (kinematic hardening)
    # Recovery    = Δt * (Rd * ( sqrt_twothirds * ϵ̇_eff ) + Rs) * κ         # recovery for kappa (isotropic hardening)
    # # σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ .* (1 - ϕ) .* Δϵ) .- ( (ϕ̇ * Δt) / (1 - ϕ) )  # deviatoric stress (trial)
    # # (^) will need this eventually | (v) this is for point simulator
    # σ̲̲⁽ᵗʳ⁾       = σ̲̲ + (2μ .* Δϵ)                                    # deviatoric stress (trial)
    # σ̲̲′⁽ᵗʳ⁾      = σ̲̲⁽ᵗʳ⁾ - (p(σ̲̲⁽ᵗʳ⁾) .* [1, 0, 0, 1, 0, 1])
    # α̲̲⁽ᵗʳ⁾       = α̲̲    .* (1 - recovery)
    # # α̲̲⁽ᵗʳ⁾       = α̲̲     * (1 - recovery)
    # κ⁽ᵗʳ⁾       = κ     * (1 - Recovery)
    # ξ̲̲⁽ᵗʳ⁾       = σ̲̲′⁽ᵗʳ⁾ - α̲̲⁽ᵗʳ⁾                      # over-stress (trial)
    # ξ_mag       = norm_symvec(ξ̲̲⁽ᵗʳ⁾)
    # # ξ_mag       = norm(ξ__)


    # # yield criterion
    # F = ξ_mag - (κ⁽ᵗʳ⁾ + β) # * (1 - ϕ)
    # # @show F
    # if F <= 0.0     # elastic
    #     # trial guesses are correct
    #     σ̲̲       = @. σ̲̲⁽ᵗʳ⁾
    #     α̲̲       = @. α̲̲⁽ᵗʳ⁾
    #     κ       =    κ⁽ᵗʳ⁾
    #     ϕ       =    ϕ
    #     ϵ̲̲      += @. Δϵ
    #     # state.ϵ_dot_plastic__    .= 0.
    # else            # plastic
    #     # Radial Return
    #     Δγ      = F / (2μ + 3\2(h + H))
    #     n̂       = ξ̲̲⁽ᵗʳ⁾ ./ ξ_mag
    #     # n̂       = ξ__ / ξ_mag
    #     σ̲̲_prev  = σ̲̲
    #     σ̲̲       = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) .* n̂
    #     α̲̲       = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) .* n̂
    #     # σ̲̲       = @. σ̲̲⁽ᵗʳ⁾ - (2μ * Δγ) * n̂
    #     # α̲̲       = @. α̲̲⁽ᵗʳ⁾ + ( h * Δγ) * n̂
    #     κ       =    κ⁽ᵗʳ⁾ + ( H * Δγ * sqrt_twothirds)
    #     # χ       = sinh(#=[=#  ( 2(2m̄ - 1) * p(σ̲̲) ) / ( (2m̄ + 1) * vonMises(σ̲̲) )  #=]=#)
    #     # ϕ       =    1 - (#={=#
    #     #         1 + (#=[=#   (1 - ϕ) ^ (1 + m̄) - 1   #=]=#) * exp(#=[=#
    #     #             (#=d̄: =# sqrt_twothirds * ϵ̇_eff) * χ * (1 + m̄) * Δt   #=]=#)
    #     #     #=}=#) ^ ( 1 / (1 + m̄) )
    #     # # (^) exact solution | (v) Forward-Euler method
    #     # # ϕ̇       = χ * (#=[=#(   1  /  ( (1 - ϕ)^m̄ )   )   -   (1 - ϕ)#=]=#) * ϵ̇_eff
    #     # # ϕ      +=    ϕ̇ * Δt
    #     # @show p(σ̲̲), vonMises(σ̲̲), χ, ϕ
    #     ϕ       = ϕ
    #     ϵ̲̲⁽ᵖ⁾   += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) ./ 2μ))
    #     # ϵ̲̲⁽ᵖ⁾   += @. (Δϵ - ((σ̲̲ - σ̲̲_prev) / 2μ))
    #     ϵ̲̲      += @. Δϵ
    # end
    # # ϵ_dot_plastic__  = @. (f * sinh(V \ (ξ_mag - κ - Y)) / ξ_mag) * ξ__
    # if ϕ > 0.99
    #     error("Element death.")
    # end
    # return σ̲̲, α̲̲, κ, ϕ, ϵ̲̲, ϵ̲̲⁽ᵖ⁾
    # # return triu_vec(σ__), triu_vec(α__), κ, triu_vec(ϵ__), triu_vec(ϵₚ__)
    # # return nothing
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
