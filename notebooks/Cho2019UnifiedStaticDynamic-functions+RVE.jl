# module B93F

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays, StructArrays
using LinearAlgebra
# using Tensors # uncomment when we can work with Tensors.jl
using DocStringExtensions

using CSV, DataFrames
using Ferrite, FerriteGmsh
using Plots, Printf
using SparseArrays, Tensors, WriteVTK
include("preamble.jl")

ContinuumMechanicsBase.I₁(x::Union{Matrix, SecondOrderTensor}) = sum(diag(x))
# ContinuumMechanicsBase.I₂(x::Union{Matrix, SecondOrderTensor}) = 2.0 \ (  ( I₁(x) ^ 2.0 )  -  ( I₁(x .^ 2.0) )  )
ContinuumMechanicsBase.I₂(x::Union{Matrix, SecondOrderTensor}) = 0.5((sum(diag(x) .^ 2.0)) + 2.0(sum(vcat(UpperTriangular(x)...) .^ 2.0)))
ContinuumMechanicsBase.I₃(x::Union{Matrix, SecondOrderTensor}) = det(x)

ContinuumMechanicsBase.I₁(x::Vector{<:Real}) = sum(x[[1, 4, 6]])
# ContinuumMechanicsBase.I₂(x::Vector{<:Real}) = 2.0 \ (  ( I₁(x) ^ 2.0 )  -  ( I₁(x .^ 2.0) )  )
ContinuumMechanicsBase.I₂(x::Vector{<:Real}) = 0.5((sum(x[[1, 4, 6]] .^ 2.0)) + 2.0(sum(x[[2, 3, 5]] .^ 2.0)))
ContinuumMechanicsBase.I₃(x::Vector{<:Real}) = det([x[1] x[2] x[3]; x[2] x[4] x[5]; x[3] x[5] x[6]])

"Maps a scalar onto the volumetric portion of the flat vector representation of a second-rank tensor."
volumetric(x::SecondOrderTensor)= x # .* diagm(diag(ones(typeof(x))))
volumetric(x::AbstractFloat)    = x .* [1, 0, 0, 1, 0, 1]
volumetric(x)                   = x .* [1, 0, 0, 1, 0, 1]

"Returns the scalar, hydrostatic portion from the flat vector representation of a second-rank tensor."
hydrostatic(x::SecondOrderTensor)= I₁(x) / 3.0
hydrostatic(x::Vector{<:Real})  = I₁(x) / 3.0
hydrostatic(x)                  = I₁(x) / 3.0

"Returns the deviatoric of the flat vector representation of a second-rank tensor."
deviatoric(x::SecondOrderTensor)= x - volumetric(SymmetricTensor{2, 3}(hydrostatic(x) .* diagm(diag(ones(typeof(x))))))
deviatoric(x::Vector{<:Real})   = x - volumetric(hydrostatic(x))
deviatoric(x)                   = x - volumetric(hydrostatic(x))

norm_symvec(tensor::Matrix) = √( sum(diag(tensor) .^ 2.0) + 2sum(vcat(UpperTriangular(tensor)...) .^ 2.0) )
norm_symvec(tensor::SecondOrderTensor) = √( sum(diag(tensor) .^ 2.0) + 2sum(vcat(UpperTriangular(tensor)...) .^ 2.0) )

"""
Structure for viscoplasticity model with loading conditions and material properties.
Here, uses the effective strain rate based on applied strain rate and loading direction.
"""
# struct Cho2019UnifiedStaticDynamic{T<:AbstractFloat} <: BammannChiesaJohnsonPlasticity.AbstractBCJMetalModel
struct Cho2019UnifiedStaticDynamicTensor{T<:AbstractFloat, S<:SymmetricTensor{4, 3, T}} <: BammannChiesaJohnsonPlasticity.AbstractBCJMetalModel
    θ       ::T         # applied temperature
    ν       ::T         # Poisson's ratio
    μ       ::T         # shear modulus
    K       ::T         # bulk modulus
    n       ::T
    ω₀      ::T
    E⁺      ::T
    V⁺      ::T
    R       ::T
    d₀      ::T
    z       ::T
    Kic     ::T
    𝒹       ::T
    𝒻       ::T
    η₀      ::T
    R₀      ::T
    P       ::T         # pressure
    ϵ̇_eff   ::T         # strain rate (effective)
    ϵₙ      ::T         # final strain
    N       ::Integer   # number of strain increments
    D⁽ᵉ⁾    ::S         # total strain tensor step
    Δt      ::T         # time step
end

"""
    $(SIGNATURES)

Outer constructor for loading conditions and material properties which assumes a Poisson's ratio of 0.5.
"""
function Cho2019UnifiedStaticDynamicTensor(Ω::BammannChiesaJohnsonPlasticity.BCJMetalStrainControl,
        ν   ::T,        # Poisson's ratio
        n   ::T,
        ω₀  ::T,
        E⁺  ::T,        # activation energy for grain growth
        V⁺  ::T,        # activation volume for grain growth
        R   ::T,        # gas constant
        d₀  ::T,        # initial grain size
        z   ::T,
        Kic ::T,        # fracture toughness
        𝒹   ::T,        # average size of second phase particles
        𝒻   ::T,        # volume fraction of second phase particles
        η₀  ::T=0.0,    # initial void nucleation density
        R₀  ::T=0.0,    # initial void radius
        P   ::T=0.0) where {T<:AbstractFloat}
    ϵ̇       = Ω.ϵ̇
    ϵₙ      = Ω.ϵₙ
    loaddir = Ω.loaddir
    M       = Ω.N + 1
    # T       = typeof(float(Ω.θ))
    # Δϵ̲̲      = zeros(T, 6)       # strain increment
    # G = E / 2(1 + ν)
    # K = E / 3(1 - 2ν)
    μ = 5.47e4    - (34.1*Ω.θ)
    K = 70000.0
    temp(i,j,k,l) = 2.0μ * (0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    D⁽ᵉ⁾ = SymmetricTensor{4, 3}(temp)
    # S       = SymmetricTensor{2, 3, T}
    # Δϵ      = zero(S) # strain increment
    Δt      = abs((ϵₙ / Ω.N) / ϵ̇)
    return Cho2019UnifiedStaticDynamicTensor(Ω.θ, ν, μ, K, n, ω₀, E⁺, V⁺, R, d₀, z, Kic, 𝒹, 𝒻, η₀, R₀, P, ϵ̇, ϵₙ, Ω.N, D⁽ᵉ⁾, Δt)
end

struct MaterialState{T, S <: SecondOrderTensor{3, T}}
    # store "converged" values
    σ̲̲   ::S # stress
    ϵ̲̲   ::S
    ϵ̲̲⁽ᵖ⁾::S # plastic strain
    α̲̲   ::S
    κ   ::T # hardening variable
    ϕ   ::T
    η   ::T
    νᵥ  ::T
    ϕ̇   ::T
    X   ::T
    XR  ::T
    XH  ::T
    Xs  ::T
    Xd  ::T
    d   ::T
end

function MaterialState(ψ)
    σ̲̲       = zero(SymmetricTensor{2, 3})       # Cauchy stress
    ϵ̲̲       = zero(SymmetricTensor{2, 3})       # total strain
    # internal state variables
    ϵ̲̲⁽ᵖ⁾    = zero(SymmetricTensor{2, 3})       # plastic strain
    α̲̲       = fill(1e-7, SymmetricTensor{2, 3}) # kinematic hardening
    κ       = 1e-7          # isotropic hardening
    ## damage
    # ϕ       = 1.0e-5       # damage
    ϕ       = 1.0e-10       # damage
    η       = ψ.η₀          # void nucleation
    # η       = p.Cnuc          # void nucleation
    # νᵥ      = 0.0           # void growth
    νᵥ      = π * (ψ.R₀^2.0)# void growth
    ϕ̇       = 0.0           # damage rate
    ## recrystallization
    X       = 1.0e-10       # total dislocation-free volume fraction
    XR      = 0.0           # total recrystallized volume fraction
    XH      = 0.0           # total reduction of recrystallized volume fraction
    Xd      = 0.5e-10       # total dynamically recrystallized volume fraction
    Xs      = 0.5e-10       # total statically recrystallized volume fraction
    d       = ψ.d₀          # average grain size
    return MaterialState(σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, ϕ, η, νᵥ, ϕ̇, X, XR, XH, Xs, Xd, d)
end

"""
Using the equations and constants from [Cho et. al. (2019)](@cite choUnifiedStaticDynamic2019), this kernel function maps the current material state and ISVs onto the next configuration.
Currently, is a literal translation of the Python code used for that publication and includes the various options for calculating recrystallization and grain growth.
Also currently includes the support for pressure-dependent systems.
"""
function update(ψ::Cho2019UnifiedStaticDynamicTensor, t, state, Δϵ̲̲, (;
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
                                    # *                         (see Eq. 4.22 in HEC dissertation).
                                    # *                         (c. f. `optimize.py` that NK=2.0 by default).
                                    # ! [20250422T1121] (JMA3): This is exponent on κ in rate equation.
                                    # !                         Bammann assumed 2 for metals for dislocation creep in Power Law
            ## torsion, tension/compression
            C₁₉, C₂₀,
            ## dynamic recrystallization
            Cx1, Cx2, Cdp,
            Cx3, Cx4, Csp,
            Cx5, Cxa, Cxb, Cxc,
            ## static RX (grain growth)
            # n, ω₀, # E⁺, V⁺, R,
            ## grain size
            # d₀, Cg1, Cg2, Cg3, z,
            Cg1, Cg2, Cg3, # z,
            # Cg1, Cg2, Cg3, Cg4, # z,
            ## damage
            ### nucleation
            # 𝒹, 𝒻, Kic, a, b, c,
            Dc, a, b, c,
            # Cnuc, Tnuc, R₀, nn, Tgrw,
            Cnuc, Tnuc, nn, Tgrw,
            ## irradiation hardening
            kr1, krt, kr2, kr3, kp1, kpt, kp2
        ); imat=0, iYS=0, tanβ₀=0.0, iREXmethod=3, iGSmethod=4, iNewton=0, kwargs...)
    # get fields from model
        σ̲̲       = state.σ̲̲
        ϵ̲̲       = state.ϵ̲̲ + Δϵ̲̲
        ϵ̲̲⁽ᵖ⁾    = state.ϵ̲̲⁽ᵖ⁾
        α̲̲       = state.α̲̲
        κ       = state.κ
        ϕ       = state.ϕ
        η       = state.η
        νᵥ      = state.νᵥ
        ϕ̇       = state.ϕ̇
        X       = state.X
        XR      = state.XR
        XH      = state.XH
        Xd      = state.Xd
        Xs      = state.Xs
        d       = state.d
        # # # @show t
        # # # @show σ̲̲
        # # # @show deviatoric(σ̲̲)
        # # # @show ϵ̲̲
        # # # @show deviatoric(ϵ̲̲)
        # # # @show ϵ̲̲⁽ᵖ⁾
        # # # @show α̲̲
        # # # @show κ
        # # # @show κₛ
        # # # @show ϕ
        # # # @show η
        # # # @show νᵥ
        # # # @show ϕ̇
    # calculation constants/functions
        sqrt_twothirds = √(2.0/3.0)
        sqrt_threehalves = √(3.0/2.0)
    # * [20250402T1345] (JMA3): Moved to `predict`
    # irradiation before damage
        Tirr = ψ.P
        M0, Si, damirr = 0.0, 0.0, 1.0
        if Tirr != 0.0
            kr = kr1 * exp(krt/Tirr)
            Si = (kr*flu) ^ (1.0/kr2)
            M0 = Kr3 * Si
            damirr = exp(  ( kp1 * exp(kpt/Tirr) * flu )  ^  ( 1.0 / kp2 )  )
        end
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
            KT0   = κ₀      + (            dK0dT * (ψ.θ-ttop)    )
            RT0   = rho0    * (  1.0  -  (   alp * (ψ.θ-ttop) )  )
            RRT0  = 1.0
            itmax = 10
            convg = 1e-12
            Niter = 0
            # Newton iterations begin
            for k in range(0, itmax)
                RRT073 = RRT0 ^ (7.0/3.0)
                RRT053 = RRT0 ^ (5.0/3.0)
                RRT023 = RRT0 ^ (2.0/3.0)
                FF = ψ.P    -    (#={=#   1.5KT0   *   (  RRT073  -  RRT053  )   *   (#=[=#
                    1.0  +  ( (3.0/4.0) * (dKdP-4.0) * (RRT023-1.0) )  #=]=#)   #=}=#)
                # define derivative of F, dF
                dF1 = ( 18.0 / 24.0 )  *  ( dKdP - 4.0 )  *  KT0  *  ( RRT0 ^ (-1.0/3.0) )  *  ( RRT073 - RRT053 )
                dF2 = 1.5KT0
                dF2 = dF2 * (  ( 7.0 / 3.0 )  *  ( RRT0 ^ (4.0/3.0) )  -  ( 5.0 / 3.0 )  *  ( RRT0 ^ (2.0/3.0) )  )
                dF2 = dF2 * (  ( 3.0 / 4.0 )  *  ( dKdP - 4.0 )  *  ( RRT023 - 1.0 )  +  1.0)
                dF  = -(dF1+dF2)
                # find corrector
                dRRT0 = -FF/dF
                # update solution
                RRT0 -= dRRT0
                # convergence check
                err = abs(dRRT0)
                err <= convg ? break : Niter += 1
                Niter >= (itmax - 1) ? println("BM convergence issue! ", err) : nothing
            end
            # pressure-temperature dependent reference density
            ρ = RRT0 * RT0
    # shear modulus
        #--- 3rd-order Finite Strain (Birch-Murnaghan EOS)
            KT0 = κ₀     +   (          dK0dT * (ψ.θ-ttop)    )
            RT0 = rho0   *   (  1.0 - (   alp * (ψ.θ-ttop) )  )
            GT0 = G0     +   (          dG0dT * (ψ.θ-ttop)    )
            b1  = (3KT0*dG0dP) - 5GT0
            b2  = 9.0(   (  KT0  ^  2.0  )   *   (  ddG0ddP  +  (
                    (1.0/KT0) * (dKdP-4.0) * dG0dP )  )   +    (  35.0GT0  /  9.0  )   )
            F   = 0.5(  ( (ρ/RT0) ^ (2.0/3.0) )  -  1.0  )
            μ   = max(  0.01,  ( (1.0+2.0F) ^ 2.5 )  *  ( GT0 + (b1*F) + 0.5b2 * (F^2.0) )  )
            ν   = 0.3
            K   = (2.0/3.0) * μ * (1.0+ν) / (1.0-2ν)
        if     imat == 1    # OFHC Cu (irradation-ISV model)
           μ = 5.47e4    - (34.1*ψ.θ)
           K = 70000.0
        elseif imat == 2    # T91 ferritic steel (Barrett et al., 2018)
           μ = 1.01e5    - (65.0*ψ.θ)
           K = 170000.0
        elseif imat == 3    # Ti6Al4V (Hukuhara&Sanpei,1993)
           μ = 4.5e4     - (20.0*ψ.θ)
           K = 85000.0
        end
        # # # @show ψ.θ, ψ.ϵₙ, ψ.N, t, μ, K
        # error("Just checking...")
    # deviatoric strain and effective strain rate
        ϵ̲̲′      = deviatoric(ϵ̲̲)
        ϵ̲̲′_mag  = norm_symvec(ϵ̲̲′)
        # ϵ̲̲′⁽ᵖ⁾   = deviatoric(ϵ̲̲⁽ᵖ⁾)
        # Δϵ̲̲⁽ᴴ⁾   = hydrostatic(ψ.Δϵ̲̲) # davg
        # Δϵ̲̲′     = deviatoric(ψ.Δϵ̲̲) # DE
        Δϵ̲̲⁽ᴴ⁾   = hydrostatic(Δϵ̲̲) # davg
        Δϵ̲̲′     = deviatoric(Δϵ̲̲) # DE
        # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
        # Δϵ̲̲̇′_mag  = norm_symvec(Δϵ̲̲′) / ψ.Δt # ddd
        # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
        # Δϵ̲̲̇′_mag  = norm_symvec(Δϵ̲̲′) / ψ.Δt # ddd
        # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
        Δϵ̲̲̇′_mag  = norm_symvec(Δϵ̲̲′) / ψ.Δt # ddd
        Δϵ̲̲̇′_mag *= sqrt_twothirds
    # trial damage
        ϕ₁⁽ᵗʳ⁾ = 1.0 - ϕ # dam1
        # @show ϕ, ϕ₁⁽ᵗʳ⁾
        ϕ₂⁽ᵗʳ⁾ = 1.0 - min(1.0, ϕ̇*ψ.Δt/ϕ₁⁽ᵗʳ⁾) # dam2
        if ϕ >= Dc
            # @info "Exceeded damage criterion. Consider this element as failed..." Dc
            σ̲̲   .= 0.0      # deviatoric stress update
            ϵ̲̲⁽ᵖ⁾ = ϵ̲̲⁽ᵖ⁾     # plastic strain update
            α̲̲   .= α̲̲        # kinematic hardening & Total strain update
            κ    = κ        # isotropic hardening update
            ϕ    = ϕ        # damage update
            η    = η
            νᵥ   = νᵥ
            ϕ̇    = ϕ̇
            X    = X        # recrystallization update
            XR   = XR
            XH   = XH
            Xd   = Xd
            Xs   = Xs
            d    = d        # grain size update
            return MaterialState(σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, ϕ, η, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d)
        end
    # hydrostatic pressure
        P = if ψ.P > 0.0
            ψ.P
        else
            # (hydrostatic(σ̲̲)*ϕ₂⁽ᵗʳ⁾) + (3.0K*Δϵ̲̲⁽ᴴ⁾*ϕ₁⁽ᵗʳ⁾)
            # hydrostatic(σ̲̲) + 3.0K*Δϵ̲̲⁽ᴴ⁾
            0.0
        end
    # deviatoric stress and invariants
        # S  = deviatoric(Sig)
        σ̲̲′  = deviatoric(σ̲̲)
        di1 = I₁(σ̲̲)
        dj2 = I₂(σ̲̲′)
        dj3 = I₃(σ̲̲′)
        JJ1 = (dj3^2.0) / (dj2^3.0)
        # # @show σ̲̲′, typeof(σ̲̲′)
        # # @show di1, dj2, dj3, JJ1
        JJ2 = (dj3    ) / (dj2^1.5)
        JJ3 = (di1    ) / (dj2^0.5)
        # # # # @show σ̲̲′
        # # # # @show di1, dj2, dj3
        # error("Just checking...")
    # temperature dependent constants
        V   = C₁ * exp(-C₂/ψ.θ)
        # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
        # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
        # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
        # Y   = C₃ * exp( C₄/ψ.θ)
        # Y   = C₃ * exp(-C₄/ψ.θ) * ((ψ.d₀/d)^ψ.z) * 0.5(1+tanh(C₁₉*(C₂₀-ψ.θ)))
        Y   = C₃ * exp( C₄/ψ.θ)
        f   = C₅ * exp(-C₆/ψ.θ)
        # # # # @show V, Y, f, C₁, C₂, C₃, C₄, C₅, C₆
        # error("Just checking...")
        # these modifiers come from Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
        if dj2 == 0.0
            djr = 1.0 - ( C₁₉ * (4.0/27.0) )
            djh = 1.0 + ( C₁₉ * (4.0/27.0) )
        else
            djr = 1.0    -    (#=[=#
                                    C₁₉   *   (  ( 4.0 / 27.0 )  -  JJ1  )
                #=]=#) - (#=[=#     C₂₀   *   (                     JJ2  )   #=]=#)
            djh = 1.0    +    (#=[=#
                                    C₁₉   *   (  ( 4.0 / 27.0 )  -  JJ1  )
                #=]=#) + (#=[=#     C₂₀   *   (                     JJ2  )   #=]=#)
        end
        # and modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
        rd  =           C₇    * exp(  -( C₈  + (1e6P*C₂₁) )  /        ψ.θ    )   *   djr
        h   = max(0.0,  C₉    * μ                                                *   djh    -    (  C₁₀  *  ψ.θ  )   )
        rs  =           C₁₁   * exp(  -( C₁₂ + (1e6P*C₂₃) )  /        ψ.θ    )
        Rd  =           C₁₃   * exp(  -( C₁₄ + (1e6P*C₂₄) )  /        ψ.θ    )   *   djr
        H   = max(0.0,  C₁₅   * μ                                                *   djh    -    (  C₁₆  *  ψ.θ  )   )
        Rs  =           C₁₇   * exp(  -( C₁₈ + (1e6P*C₂₆) )  /        ψ.θ    )
        # # # # @show C₁₅, μ, djh
        # # # # @show C₁₅ * μ * djh
        # # # # @show djr, djh, rd, h, rs, Rd, H, Rs, Rdc
        # error("Just checking...")
    # yield surface parameters
        # iYS: 0-Pressure insensitive (Mises);
        #      1-Pressure sensitive (Shear-Mises);
        #      2-Pressure sensitive (TANH)
        if     iYS == 0
            Yₚ = 0.
        elseif iYS == 1
            tanB = tanβ₀
            Pa = Pₖ₁   *   (  ( 1.0 + exp(-Pₖ₂/ψ.θ) )  ^  ( -Pₖ₃ )  )
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
            #Yp = Pk1*exp(-Pk2*ψ.θ)*tanh(B2*P[i])
            #Yp = (Pk1*(1. + exp(-Pk2/ψ.θ))^(-Pk3))*tanh(B2*P[i])
            β₁ = max(  1e-10,    ( Pₖ₁ - (Pₖ₂*ψ.θ) )  )
            Yₚ = β₁  *  tanh( (Pₖ₃/β₁) * P )
            # # # # @show β₁, ψ.θ, Pₖ₁, Pₖ₂, Pₖ₃, P, Yₚ
            # error("Just checking...")
        else
            error("iYS > 2 which is not supported.")
        end
    # viscous stress
        # β = V      *      log(     f     \     (#={=#
        #             Δϵ̲̲̇′_mag                               +    sqrt(#=[=#
        #             (    Δϵ̲̲̇′_mag               ^  2.0  )  +  (  f  ^  2.0  )   #=]=#)
        #     #=}=#)     )
        #... Using sinh^n for strain rate-stress curve's smooth connection
        # β = V      *      log(     f     \     (#={=#
        #         (   Δϵ̲̲̇′_mag   ^   (  1.0  /  NK  )   )    +    sqrt(#=[=#
        #             (  ( Δϵ̲̲̇′_mag ^ (1.0/NK) )  ^  2.0  )  +  (  f  ^  2.0  )   #=]=#)
        #     #=}=#)     )
        # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
        # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
        # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
        # β = V * asinh(sqrt_twothirds*Δϵ̲̲̇′_mag/f)
        # β = V * asinh(Δϵ̲̲̇′_mag/f)
        β = V * asinh(Δϵ̲̲̇′_mag/f)
        # # # # @show Yₚ, Be, V, Δϵ̲̲̇′_mag, NK, f
        # error("Just checking...")
    # previous alpha magnitude
        # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
        # α̲̲_mag   = t <= 1.01ψ.Δt ? 0.0 : norm_symvec(α̲̲)
        # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
        # α̲̲_mag   = t <= 1.01ψ.Δt ? 0.0 : norm_symvec(α̲̲)
        # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
        α̲̲_mag   = t <= 1.01ψ.Δt ? 0.0 : norm_symvec(α̲̲)
        α̲̲_mag  *= sqrt_threehalves
        # α̲̲_mag  /= sqrt_threehalves
        # # # @show α̲̲
        # # # @show sqrt_threehalves, α̲̲_mag
        # # α̲̲_mag  /= sqrt_threehalves
        # # # @show α̲̲_mag
        # # error("Just checking...")
        # if t > ψ.Δt
        #     error("Just checking...")
        # end
    # REX Model
        ## REX calculation: separated DRX and SRX equations
            if     iREXmethod == 0
                xx  = 0.0
                dXR = 0.0
                dXH = 0.0
                dXd = 0.0
                dXs = 0.0
            elseif iREXmethod == 1 # Euler Method (explicit)
                # KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                # dAlpha  = (  h  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                # dAlpha  = max(0.0, dAlpha)
                # dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^ 2.0 )  )
                # dKappa  = max(0.0, dKappa)
                # KAlMu1  = μ  \  ( dKappa + dAlpha )
                # Cxd     = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  *  ψ.Δt  )
                # Cxs     = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  ψ.θ  )   *   (  KAlMu             *  ψ.Δt  )
                # Ch      = Cx5 * KAlMu1 * ψ.Δt
                KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                dAlpha  = (  h  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^ 2.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMuX  = κ  \  ( dKappa + dAlpha )
                Cxd     = Cx1   *   exp(  -( Cx2 + (P*1e6Cdp) )  /  ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  *  ψ.Δt  )
                Cxs     = Cx3   *   exp(  -( Cx4 + (P*1e6Csp) )  /  ψ.θ  )   *   (  KAlMu             *  ψ.Δt  )
                Ch      = Cx5 * KAlMuX * ψ.Δt
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
                # # ? [20250402T1149] (JMA3): Maybe this (v) should be included? It's not originally...
                # # * [20250402T1151] (JMA3): Maybe this section of reassignment is redundant
                # # * [20250424T0955] (JMA3): I commented this out because it seems redundant.
                # # * ========================================================================
                # # *                         since these get updated at the end anyway.
                # # XR     += dXR # ! update ISV
                # XH     += dXH # ! update ISV
                # Xd     += dXd # ! update ISV
                # Xs     += dXs # ! update ISV
                # # * ========================================================================
                # # dX      = dXR - dXH
                xx      = X + dX
            elseif iREXmethod == 2 # explicit exponential integration algorithm
                KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                # dAlpha  = (  h  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                # dAlpha  = max(0.0, dAlpha)
                # dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^ 2.0 )  )
                # dKappa  = max(0.0, dKappa)
                # KAlMu1  = μ  \  ( dKappa + dAlpha )
                # Cxd     = Cx1   *   exp(  -( Cx2 + (   P*Cdp) )  /        ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  *  ψ.Δt  )
                # Cxs     = Cx3   *   exp(  -( Cx4 + (   P*Csp) )  /        ψ.θ  )   *   (  KAlMu             *  ψ.Δt  )
                # Ch      = Cx5 * KAlMu1 * ψ.Δt
                KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                dAlpha  = (  h  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                dAlpha  = max(0.0, dAlpha)
                dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^ 2.0 )  )
                dKappa  = max(0.0, dKappa)
                KAlMuX  = κ  \  ( dKappa + dAlpha )
                Cxd     = Cx1   *   exp(  -( Cx2 + (   P*1e6Cdp) )  /        ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  *  ψ.Δt  )
                Cxs     = Cx3   *   exp(  -( Cx4 + (   P*1e6Csp) )  /        ψ.θ  )   *   (  KAlMu             *  ψ.Δt  )
                Ch      = Cx5 * KAlMuX * ψ.Δt
                Udt     = ( Cxd + Cxs )  *  ( X ^ Cxa )  *  ( (1.0-X) ^ (Cxb-1.0) )
                Udt     = Udt   +   (  Ch  *  (  X ^ (Cxc-1.0)  )  )
                Vdt     = ( Cxd + Cxs )  *  ( X ^ Cxa )  *  ( (1.0-X) ^ (Cxb-1.0) )
                xx      = (  X   *   exp( -Udt )  )   +   (  Vdt  *  ( (1.0-exp(-Udt)) / Udt )  )
                pX0     = Cxd + Cxs
                dXR     = pX0          *  ( X ^ Cxa )  *  ( (1.0-X) ^  Cxb  )
                dXH     = Ch * (X^Cxc)
            elseif iREXmethod == 3 # RK4-explicit method
                # KAlMu   = μ  \  ( (κ^1.0) + ((                 α̲̲_mag)^1.0) )
                KAlMu   = μ  \  ( (κ^2.0) + ((sqrt_threehalves*α̲̲_mag)^2.0) )
                # KAlMu   = μ  \  ( (κ^2.0) + ((                 α̲̲_mag)^2.0) )
                # # KAlMu   = μ  \  ( ((κ/ψ.Δt)^2.0) + (((sqrt_threehalves*α̲̲_mag)/ψ.Δt)^2.0) )
                # KAlMu   = μ  \  ( κ + (sqrt_threehalves*α̲̲_mag) )
                # # KAlMu   = μ  \  ( (κ^2.0) + ((sqrt_threehalves*α̲̲_mag)^2.0) )
                # dAlpha  = (  h  *  (               Δϵ̲̲̇′_mag)  )   -   (  ( (sqrt_twothirds* rd*(               Δϵ̲̲̇′_mag)) + rs )  *  ( (                 α̲̲_mag) ^ 3.0 )  )
                # # dAlpha  = (  h  *  (sqrt_twothirds*Δϵ̲̲̇′_mag)  )   -   (  ( (sqrt_twothirds* rd*(sqrt_twothirds*Δϵ̲̲̇′_mag)) + rs )  *  ( (sqrt_threehalves*α̲̲_mag) ^  NK )  )
                # # dAlpha  = (  h  *  ( 1 - X )  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^  NK )  )
                # # dAlpha  = (  h  *  ( 1 - X )  *  (Δϵ̲̲̇′_mag*ψ.Δt)  )   -   (  ( (sqrt_twothirds* rd*(Δϵ̲̲̇′_mag*ψ.Δt)) + rs )  *  ( α̲̲_mag ^  NK )  )
                # # dAlpha  = (  h  *  Δϵ̲̲̇′_mag  )   -   (  ( (                rd*Δϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^  NK )  )
                dAl = [(h*Δϵ̲̲′[k]/ψ.Δt - ((rd*Δϵ̲̲̇′_mag+rs)*α̲̲_mag*α̲̲[k])) for k in [1, 4, 6, 2, 3, 5]]
                # # # # @show h, ψ.Δt, rd, Δϵ̲̲̇′_mag, rs, α̲̲_mag
                # # # # @show Δϵ̲̲′
                # # # # @show α̲̲
                # error("Just checking...", dAl)
                dAlpha  = sum(dAl[[1, 4, 6]] .^ 2.0)
                dAlpha += sum(2.0 .* (dAl[[2, 3, 5]] .^ 2.0))
                dAlpha  = √(3dAlpha/2)
                dAlpha  = max(0.0, dAlpha) # * ((ψ.d₀/d)^ψ.z)
                # dKappa  = (  H  *  (               Δϵ̲̲̇′_mag)  )   -   (  ( (sqrt_twothirds* Rd*(               Δϵ̲̲̇′_mag)) + Rs )  *  (     κ ^ 3.0 )  )
                dKappa  = (  H  *  (               Δϵ̲̲̇′_mag)  )   -   (  ( (                Rd*(               Δϵ̲̲̇′_mag)) + Rs )  *  (     κ ^  NK )  )
                # dKappa  = (  H  *  (sqrt_twothirds*Δϵ̲̲̇′_mag)  )   -   (  ( (sqrt_twothirds* Rd*(sqrt_twothirds*Δϵ̲̲̇′_mag)) + Rs )  *  (     κ ^  NK )  )
                # dKappa  = (  H  *  ( 1 - X )  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^  NK )  )
                # # dKappa  = (  H  *  ( 1 - X )  *  (Δϵ̲̲̇′_mag*ψ.Δt)  )   -   (  ( (sqrt_twothirds* Rd*(Δϵ̲̲̇′_mag*ψ.Δt)) + Rs )  *  (     κ ^  NK )  )
                # # dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (                Rd*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^  NK )  )
                dKappa  = max(0.0, dKappa) # * ((ψ.d₀/d)^ψ.z)
                KAlMuX  = μ  \  ( dKappa + dAlpha )
                # KAlMuX  = μ  \  ( ((dKappa)^2.0) + ((dAlpha)^2.0) )
                # Cxd     = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /  ( ψ.R * ψ.θ )  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  )
                Cxd     = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /        ψ.θ    )   *   (  KAlMu  *  (               Δϵ̲̲̇′_mag)  )
                # Cxd     = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /        ψ.θ    )   *   (  KAlMu  *  (sqrt_twothirds*Δϵ̲̲̇′_mag)  )
                Cxs     = Cx3   *   exp(  -( Cx4 + (   P*Csp) )  /        ψ.θ    )   *      KAlMu
                Ch      = Cx5 * KAlMuX
                # Cxd     = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /        ψ.θ    )   *   (  KAlMuX  *  Δϵ̲̲̇′_mag  )
                # Cxs     = Cx3   *   exp(  -( Cx4 + (   P*Csp) )  /        ψ.θ    )   *      KAlMuX
                # Ch      = Cx5 * KAlMu
                CC      = Cxd + Cxs
                k₁      = CC    *    (                 X                          ^   Cxa       )    *    (   (  1.0  -    X                 )   ^   Cxb   )
                k₁      = k₁    -    (   Ch   *   (    X                          ^   Cxc   )   )
                k₂      = CC    *    (            (    X  + 0.5( ψ.Δt * k₁ )  )   ^   Cxa       )    *    (   (  1.0  -  ( X + 0.5(ψ.Δt*k₁) )  )   ^   Cxb   )
                k₂      = k₂    -    (   Ch   *   (  ( X  + 0.5( ψ.Δt * k₁ )  )   ^   Cxc   )   )
                # # # # @show i, d
                # # # # @show μ, κ, sqrt_threehalves*α̲̲_mag
                # # # # @show KAlMu, dAl
                # # # # @show CC, X, ψ.Δt, k₂, Cxa, Cxb
                # # # # @show (X + 0.5*ψ.Δt*k₂)
                # # # # @show (X + 0.5*ψ.Δt*k₂)^Cxa
                # # # # @show 1. - (X + 0.5*ψ.Δt*k₂)
                # # # # @show CC*((X + 0.5*ψ.Δt*k₂)^Cxa)*(1. - (X + 0.5*ψ.Δt*k₂))^Cxb
                # error("Just checking...")
                k₃      = CC    *    (            (    X  + 0.5( ψ.Δt * k₂ )  )   ^   Cxa       )    *    (   (  1.0  -  ( X + 0.5(ψ.Δt*k₂) )  )   ^   Cxb   )
                k₃      = k₃    -    (   Ch   *   (  ( X  + 0.5( ψ.Δt * k₂ )  )   ^   Cxc   )   )
                k₄      = CC    *    (            (    X  +    ( ψ.Δt * k₃ )  )   ^   Cxa       )    *    (   (  1.0  -  ( X +    (ψ.Δt*k₃) )  )   ^   Cxb   )
                k₄      = k₄    -    (   Ch   *   (  ( X  +    ( ψ.Δt * k₃ )  )   ^   Cxc   )   )
                xx      = X   +   (  1.0  /  6.0  )   *   (  ( k₁ + 2.0(k₂+k₃) + k₄ )  *  ψ.Δt  )
                # Cxd    *= ψ.Δt
                # Cxs    *= ψ.Δt
                # Ch     *= ψ.Δt
                pX0     = Cxd  +  Cxs
                dXd     = Cxd  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXs     = Cxs  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXR     = pX0  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                dXH     =  Ch  *  ( xx ^ Cxc )
                # # # # @show KAlMu, dKappa, dAlpha, KAlMu1
                # # # # @show Cxd, Cxs, Ch, pX0
                # # # # @show dXd, dXs, dXR, dXH
                # error("Just checking...")
            elseif iREXmethod >= 4 # implicitly solve functions using Newton-Rapson method
                Nitmax = 20
                Ntol   = 1e-6
                xx     = 0.5
                for k in range(0, Nitmax)
                    if     iREXmethod == 4 # Euler Method (implicit)
                        KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                        dAlpha  = (  h  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) + rs )  *   ( α̲̲_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        # dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) + Rs )  *   (     κ ^ 2.0 )  )
                        dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) + Rs )  *   (     κ ^  NK )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMuX  = μ  \  ( dKappa + dAlpha )
                        Cxd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  *  ψ.Δt  )
                        Cxs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  ψ.θ  )   *   (  KAlMu             *  ψ.Δt  )
                        Ch      = Cx5 * KAlMuX * ψ.Δt
                        # TODO [20250331T1119] (JMA3): come back to decrement this section instead
                        F       = X  +   (  ( Cxd + Cxs )  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )  )
                        F       = F  - ( Ch * (xx^Cxc) ) - xx
                        dF      =    - ( Cxd + Cxs )  *  Cxb  *  ( (1.0-xx) ^ (Cxb-1.0) )  *  ( xx ^  Cxa )
                        dF      = dF + ( Cxd + Cxs )  *  Cxa  *  ( (1.0-xx) ^  Cxb      )  *  ( xx ^ (Cxa-1.0))
                        dF      = dF - (      Ch    *  Cxc  *  ( xx ^ (Cxc-1.0) )  )
                        dF      = dF - 1.0
                    elseif iREXmethod == 5 # exponential integration algorithm (asymptotic)
                        KAlMu   = μ  \  ( (κ^2.0) + (α̲̲_mag^2.0) )
                        dAlpha  = (  h  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) + rs )  *  ( α̲̲_mag ^ 2.0 )  )
                        dAlpha  = max(0.0, dAlpha)
                        # dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^ 2.0 )  )
                        dKappa  = (  H  *  Δϵ̲̲̇′_mag  )   -   (  ( (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) + Rs )  *  (     κ ^  NK )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMuX  = μ  \  ( dKappa + dAlpha )
                        Cxd      = Cx1   *   exp(  -( Cx2 + (P*Cdp) )  /  ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  *  ψ.Δt  )
                        Cxs      = Cx3   *   exp(  -( Cx4 + (P*Csp) )  /  ψ.θ  )   *   (  KAlMu             *  ψ.Δt  )
                        Ch      = Cx5 * KAlMuX * ψ.Δt
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
                    elseif iREXmethod == 6 # exponential integration algorithm (trapezoidal)
                        KAlMu   = μ  \  ( (κ^2.0) + (                  α̲̲_mag ^2.0) )
                        # KAlMu   = μ  \  ( (κ^2.0) + ((sqrt_threehalves*α̲̲_mag)^2.0) )
                        dAlpha  = (  h  *                  Δϵ̲̲̇′_mag   )   -   (  ( (sqrt_twothirds* rd*(               Δϵ̲̲̇′_mag)) + rs )  *  ( (                 α̲̲_mag) ^ 2.0 )  )
                        # dAlpha  = (  h  *  (sqrt_twothirds*Δϵ̲̲̇′_mag)  )   -   (  ( (sqrt_twothirds* rd*(sqrt_twothirds*Δϵ̲̲̇′_mag)) + rs )  *  ( (sqrt_threehalves*α̲̲_mag) ^  NK )  )
                        dAlpha  = max(0.0, dAlpha)
                        dKappa  = (  H  *  (               Δϵ̲̲̇′_mag)  )   -   (  ( (sqrt_twothirds* Rd*(               Δϵ̲̲̇′_mag)) + Rs )  *  (     κ ^ 2.0 )  )
                        # dKappa  = (  H  *  (sqrt_twothirds*Δϵ̲̲̇′_mag)  )   -   (  ( (sqrt_twothirds* Rd*(sqrt_twothirds*Δϵ̲̲̇′_mag)) + Rs )  *  (     κ ^  NK )  )
                        dKappa  = max(0.0, dKappa)
                        KAlMuX  = μ  \  ( dKappa + dAlpha )
                        Cxd     = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /  ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  *  ψ.Δt  )
                        Cxs     = Cx3   *   exp(  -( Cx4 + (1e6P*Csp) )  /  ψ.θ  )   *   (  KAlMu              *  ψ.Δt  )
                        Ch      = Cx5 * KAlMuX * ψ.Δt
                        # Cxd     = Cx1   *   exp(  -( Cx2 + (1e6P*Cdp) )  /  ψ.θ  )   *   (  KAlMu  *  Δϵ̲̲̇′_mag  ) # *  ψ.Δt  )
                        # Cxs     = Cx3   *   exp(  -( Cx4 + (1e6P*Csp) )  /  ψ.θ  )   *   (  KAlMu              ) # *  ψ.Δt  )
                        # Ch      = Cx5 * KAlMuX # * ψ.Δt
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
                    # dxx/= ψ.Δt
                    xx  = max(1e-6, min(0.9999999, xx + dxx))
                    # xx  = max(1e-6, min(0.9999999, xx + dxx*ψ.Δt))

                    pX0     = Cxd  +  Cxs
                    dXd     = Cxd  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                    # dXd    /= ψ.Δt
                    dXs     = Cxs  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                    # dXs    /= ψ.Δt
                    dXR     = pX0  *  ( xx ^ Cxa )  *  ( (1.0-xx) ^ Cxb )
                    # dXR    /= ψ.Δt
                    dXH     =  Ch  *  ( xx ^ Cxc )
                    # dXH    /= ψ.Δt

                    if abs(dxx) <= Ntol
                    # if abs(dxx*ψ.Δt) <= Ntol
                        break
                    end
                    if k >= Nitmax-1
                        println("Newton-Rapson Convergence Issue: k >= Nitmax")
                    end
                end
            else
                error("iREXMethod > 6 which is not supported.")
            end

            # Final solution and estimate volume of DRX and SRX
            xx  = max(1e-10, min(0.9999999, xx))
            dX  = xx - X
            X0  = 1.0  -  ( dX / (1.0-X) )
            X   = xx # ! update ISV
            (dX < 0.0)  ?  (X0 = 1.0)  :  nothing # ∵ REX is irreversible
            XR += dXR * ψ.Δt # ! update ISV
            XH += dXH * ψ.Δt # ! update ISV
            #dXd   = dXR*(Cd/pX0)
            #dXs   = dXR*(Cs/pX0)
            Xd += dXd * ψ.Δt # ! update ISV
            Xs += dXs * ψ.Δt # ! update ISV
            Rx  = (1.0-X) ^ NK
            # # # # @show xx, dX, X0, X
            # # # # @show XR, XH, Xd, Xs
            # error("Just checking...")
        ## Grain size kinetics (SGG and grain refinement rate)
            # Grain size rate integration method
            # 0-explicit; 1-implicit; 2-analytic; 3-earlier model (IJP,2019)
            # iGSmethod = 0 # [20250402T1523] (JMA3): I commented this out to let the keyword argument have precedence.
            dim1 = d
            if     iGSmethod == 0
                d = dim1
            elseif iGSmethod == 1 # Forward Euler (explicit)
                # static grain growth rate
                dg      = dim1
                ω    = ψ.ω₀   *   exp(  -( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )
                dsgg    = ω   /   (  ψ.n  *  ( dg ^ (ψ.n-1.0) )  )
                # dynamic grain size reduction rate (new version: EPSL2020)
                dred    = Cg1 * X * Δϵ̲̲̇′_mag * (dg^Cg2)
                # total grain size change rate
                d       = dg  +  ( (dsgg-dred) * ψ.Δt ) # ! update ISV
                # Z       = ddd*exp((sxE + P[i]*1.e6*sxV)/(ψ.R*ψ.θ))
                # dss     = (sxk/(Cg3*sxn*0.3))^(1./(sxn-1.+Cg2))*Z^(-(1./(sxn-1.+Cg2)))
                # # # # @show dr, dsgk, dsgg, dred, d
                # error("Just checking...")
            elseif iGSmethod == 2 # Explicit RK4
                dg      = d
                ω    = ψ.ω₀  *  exp( -ψ.E⁺ / (ψ.R*ψ.θ) )
                Grd0    = Cg1 * Xd * Δϵ̲̲̇′_mag
                gk1     = ω   /   (  ψ.n  *  (  dg            ^ (ψ.n-1.0) )  )   -   (  Grd0  *  (  dg            ^ Cg2 )  )
                gk2     = ω   /   (  ψ.n  *  ( (dg+0.5gk1*ψ.Δt) ^ (ψ.n-1.0) )  )   -   (  Grd0  *  ( (dg+0.5gk1*ψ.Δt) ^ Cg2 )  )
                gk3     = ω   /   (  ψ.n  *  ( (dg+0.5gk2*ψ.Δt) ^ (ψ.n-1.0) )  )   -   (  Grd0  *  ( (dg+0.5gk2*ψ.Δt) ^ Cg2 )  )
                gk4     = ω   /   (  ψ.n  *  ( (dg+   gk3*ψ.Δt) ^ (ψ.n-1.0) )  )   -   (  Grd0  *  ( (dg+   gk3*ψ.Δt) ^ Cg2 )  )
                d       = dg  +  ( 1.0 / 6.0 )  *  ( gk1 + 2.0(gk2+gk3) + gk4 )  *  ψ.Δt
            elseif iGSmethod == 3 # Backward Euler: a = 1 (implicit); a = 0.5 (Crank-Nicholson)
                λ       = 0.5
                Nitmax  = 20
                Convg   = 1e-6
                dg      = dim1
                # dsgk    = sxk*exp(-(sxE + P[i]*1.e6*sxV)/(ψ.R*ψ.θ))
                # time downscaling factor for matching to n=4
                # tscl    = t  ^  ( (ψ.n/4.0) - 1.0 )
                # ω       = ψ.ω₀   *   exp(  -( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )   *   tscl
                # ω       = ψ.ω₀   *   exp(  -( ψ.E⁺             )  /  ( ψ.R * ψ.θ )  )
                ω       = ψ.ω₀   *   exp(  -( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )
                xx      = dg
                for k in range(1, Nitmax)
                    F   = dg      +      (#={=#     ω     *     ψ.Δt     /     (#=[=#    ψ.n    *    (
                            (  ( 1.0 - λ )  *  ( dg ^ (ψ.n-1.0) )  )   +   (  λ  *  ( xx ^ (ψ.n-1.0) )  )
                        )    #=]=#)     #=}=#)
                    # # # @show dg, ω, ψ.Δt, ψ.n, λ, xx, Cg1, X, Δϵ̲̲̇′_mag, Cg2
                    # # # @show Cg1 * X * Δϵ̲̲̇′_mag * ψ.Δt
                    # # # @show dg^Cg2
                    # # # @show ( (1.0-λ) * (dg^Cg2) )
                    # # # @show xx^Cg2
                    # # # @show ( λ * (xx^Cg2) )
                    # # # @show ( (1.0-λ) * (dg^Cg2) ) + ( λ * ((xx+0im)^Cg2) )
                    # # # @show (Cg1 * X * Δϵ̲̲̇′_mag * ψ.Δt) * (( (1.0-λ) * (dg^Cg2) ) + ( λ * ((xx+0im)^Cg2) ))
                    F  -= Cg1   *   X   *   Δϵ̲̲̇′_mag   *   ψ.Δt   *   (
                            ( (1.0-λ) * (dg^Cg2) )  +  ( λ * (xx^Cg2) )  )
                    # F  -= Cg1   *   X   *   (sqrt_twothirds*Δϵ̲̲̇′_mag)   *   ψ.Δt   *   (
                    #         ( (1.0-λ) * (dg^Cg2) )  +  ( λ * (xx^Cg2) )  )
                    F  -= xx
                    dF  = (  ( ω * ψ.Δt * λ * (1.0-ψ.n) / ψ.n )  *  ( xx ^ -ψ.n )  )   -   1.0
                    dF -= Cg1    *    X    *    Δϵ̲̲̇′_mag    *    ψ.Δt    *    (
                            (  Cg2 * λ * ( (xx+0im) ^ (Cg2-1.0) )  )   )
                    # dF -= Cg1    *    X    *    (sqrt_twothirds*Δϵ̲̲̇′_mag)    *    ψ.Δt    *    (
                    #         (  Cg2 * λ * ( (xx+0im) ^ (Cg2-1.0) )  )   )
                    dxx = -F / dF
                    xx += dxx

                    if abs(dxx) <= Convg
                    # if abs(dxx/xx) <= Convg
                        break
                    end
                    if k >= Nitmax
                        println("N-R Convg Issue for Grain Size: k >= Nitmax", dxx, k)
                    end
                end
                d       = abs(xx) # ! update ISV
                # # prefct  = ( ψ.ω₀ * tscl / (Cg1*ψ.n*X) )  ^  ( 1.0 / (ψ.n-1.0+Cg2) )
                # # dsss    = prefct     *     (    (#=[=#
                # #         Δϵ̲̲̇′_mag   *   exp(  ( ψ.E⁺ + (1e6P*ψ.V⁺) ) / ( ψ.R * ψ.θ )  )
                # #     #=]=#)    ^    (   -1.0   /   (  ψ.n  -  1.0  +  Cg2  )   )    )
                # prefct  = ( ψ.ω₀ / (Cg1*ψ.n*X) )  ^  ( 1.0 / (ψ.n-1.0+Cg2) )
                # dsss    = prefct     *     (    (#=[=#
                #         Δϵ̲̲̇′_mag   *   exp(  ( ψ.E⁺             ) / ( ψ.R * ψ.θ )  )
                #     #=]=#)    ^    (   -1.0   /   (  ψ.n  -  1.0  +  Cg2  )   )    )
            elseif iGSmethod == 4 # analytical solution
                # static grain growth
                ω    = ψ.ω₀   *   exp(  -( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )
                # ! update ISV
                # // ? [20250401T1206] (JMA3): what is `d0`
                # [20250422T1126] (JMA3): `d0` is the initial grain size.
                d       = ψ.d₀    +    (   ω   *   t   *   (  t  ^  ( (ψ.n/4.0) - 1.0 )  )   )    ^    (   1.0   /   ψ.n   )
            elseif iGSmethod == 5 # Initial model (Cho et al. (2019) IJP)
                dg  = d
                # time downscaling factor for matching to n=4
                # tscl= t  ^  ( (ψ.n/4.0) - 1.0 )
                # ω   = ψ.ω₀   *   exp(  -( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )   *   tscl
                ω   = ψ.ω₀   *   exp(  -( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )
                # ω   = ψ.ω₀   *   exp(  -( ψ.E⁺               )  /  ( ψ.R * ψ.θ )  )
                # ? [20250331T1347] (JMA3): What even is this `if`-statement?
                # * [20250729T1140] (JMA3): I guess it sets the steady-state grain size for creep
                dss = if Δϵ̲̲̇′_mag == 0.0
                    dg
                else
                    Z =                 Δϵ̲̲̇′_mag    *   exp(  ( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )
                    # Z =                 Δϵ̲̲̇′_mag    *   exp(  ( ψ.E⁺               )  /  ( ψ.R * ψ.θ )  )
                    # Z = (sqrt_twothirds*Δϵ̲̲̇′_mag)   *   exp(  ( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )
                    # # @show Z, Cg1, Cg2, Cg1 * (Z^-Cg2)
                    # error("Just checking...")
                    Cg1 * (Z^-Cg2)
                end
                # ? [20250331T1350] (JMA3): Why the addition, subtraction, and increment?
                # ddgrw   = (  ( (ω*ψ.Δt) + (dg^ψ.n) )  ^  ( 1.0 / ψ.n )  )   -   dg
                # dg     += ddgrw
                # # @show ω, ψ.n, dg,  ω   /   (  ψ.n  *  ( dg ^ (ψ.n-1.0) )  )
                # error("Just checking...")
                ddgrw   = ω   /   (  ψ.n  *  ( dg ^ (ψ.n-1.0) )  )
                dg     += ddgrw*ψ.Δt
                ds      = min(dss, dg)
                # # @show ddgrw, dg, ds
                # error("Just checking...")
                # ddred   = -Cg3  *  Xd  *  Δϵ̲̲̇′_mag  *  ψ.Δt  *  dg  *    (dg-ds)
                # ddred   = Cg3  *  dXd  *  Δϵ̲̲̇′_mag  *  ψ.Δt  *  dg  *  ( (ds-dg) ^ 2.0 )
                ddred   = Cg3  *  dXd  *  dg  *  ( (ds-dg) ^ 2.0 )
                # ddred   = Cg3  *  (Xd/ψ.Δt)  *  dg  *  ( (ds-dg) ^ 2.0 )  *  ψ.Δt
                # ddred   = max((ds-dg), ddred)
                ddred   = max((ds-dg), ddred*ψ.Δt)
                d       = dg - ddred*ψ.Δt # ! update ISV
                # ds     = min(dss, dg)
                # # ω   = ψ.ω₀   *   exp(  -( ψ.E⁺ + (1e6P*ψ.V⁺) )  /  ( ψ.R * ψ.θ )  )
                # ω   = ψ.ω₀   *   exp(  -( ψ.E⁺               )  /  ( ψ.R * ψ.θ )  )
                # # Grd0= Cg3 * (Xd^Cg4) * Δϵ̲̲̇′_mag
                # Grd0= Cg3 * dXd # * Δϵ̲̲̇′_mag
                # gk1 = ω   /   (  n  *  (   dg ^ (n-1.0) )  )   -   (  Grd0  *    dg  *  ( (ds-dg) ^ 2.0 )  )
                # # # @show dsgk, n, dg, Grd0, ds, gk1
                # # # error("Just checking...")
                # # if t > ψ.Δt
                # #     error("Just checking...")
                # # end
                # gkgt= dg + (0.5gk1*ψ.Δt)
                # gk2 = ω   /   (  n  *  ( gkgt ^ (n-1.0) )  )   -   (  Grd0  *  gkgt  *  ( (ds-gkgt) ^ 2.0 )  )
                # gkgt= dg + (0.5gk2*ψ.Δt)
                # gk3 = ω   /   (  n  *  ( gkgt ^ (n-1.0) )  )   -   (  Grd0  *  gkgt  *  ( (ds-gkgt) ^ 2.0 )  )
                # gkgt= dg + (   gk3*ψ.Δt)
                # gk4 = ω   /   (  n  *  ( gkgt ^ (n-1.0) )  )   -   (  Grd0  *  gkgt  *  ( (ds-gkgt) ^ 2.0 )  )
                # d   = dg  +  ( 1.0 / 6.0 )  *  ( gk1 + 2.0(gk2+gk3) + gk4 )  *  ψ.Δt # ! update ISV
                # # # @show gk1, gk2, gk3, gk4
                # # # @show dg, dsgk, Grd0, d
                # # # error("Just checking...")
                # # if t > ψ.Δt
                # #     error("Just checking...")
                # # end
            else
                error("iGSmethod > 4 not supported")
            end
        ## Hall-Petch effect
            idzz = 0
            # // ? [20250401T1206] (JMA3): what is `d0`
            # [20250422T1126] (JMA3): `d0` is the initial grain size.
            dzz1, dzz0 = if idzz == 0
                # # @show ψ.d₀, d, ψ.d₀/d
                # error("Just checking...")
                ( (ψ.d₀/d) ^ ψ.z,             1.0 )
            elseif idzz == 1
                (         1.0,     (dim1/d) ^ ψ.z )
            elseif idzz == 2
                ( (ψ.d₀/d)     ,   (dim1/d)     ) .^ ψ.z
            else
                error("idzz > 2 which is not supported.")
            end
            # # # # @show idzz, dzz1, dzz0
            # error("Just checking...")
            # [20250401T1042] (JMA3): these comments (v) are from HEC's original code
            # d0 = 1. !Turn on if absolute grain size-stress relation is used
            # YT  = YT*dzz1 ! Turn on if grain size dependent yield is used
    # elastic prediction
        twoμ = 2.0μ
        #--- trial deviatoric stress
            # σ̲̲′⁽ᵗʳ⁾ = ( σ̲̲′ .* ϕ₂⁽ᵗʳ⁾ )  +  ( twoμ .* Δϵ̲̲′ .* ϕ₁⁽ᵗʳ⁾) # Str
            σ̲̲′⁽ᵗʳ⁾ = ( σ̲̲′ .* ϕ₂⁽ᵗʳ⁾ )  +  ( deviatoric(ψ.D⁽ᵉ⁾ ⊡ (ϵ̲̲ - ϵ̲̲⁽ᵖ⁾)) .* ϕ₁⁽ᵗʳ⁾) # Str
        #--- use of Newton Method for DG and Kappa
        # iNewton = 0 # [20250422T1128] (JMA3): Commented out to let keyword argument have precendence.
        # #--- irradiation hardening effect
            # Hir = (1.0+Si) ^ 2.0
        #--- trial kappa
            # if iNewton == 0
            # end
            rdrsk   = 1.0   +   (  ( Rs + (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  (     κ ^ (NK-1.0) )  *  dzz1  )
            κ⁽ᵗʳ⁾   = κ * (X0*dzz0/rdrsk) # Ktr
            # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
            # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
            # # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
            # # rdrsk   = 1.0   -   (  ( Rs + (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  (     κ ^ (NK-1.0) )  )
            # # rdrsk   = 1.0   -   (  ( Rs + (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  (     κ ^ (NK-1.0) )  *  dzz1  )
            # rdrsk   = 1.0   -   (  ( Rs + (sqrt_twothirds* Rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  (     κ ^ (NK-1.0) )  *  dzz1  )
            # κ⁽ᵗʳ⁾   = κ * (X0*dzz0*rdrsk) # Ktr
        # #--- trial M in isotropic hardening (output only)
            # rdrssk  = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  ( κₛ ^ (NK-1.0) ) *  dzz1  )
            # # Kstr    = κₛ * X0 * dzz0 / rdrssk
            # κₛ⁽ᵗʳ⁾  = κₛ * X0 * dzz0 / rdrssk
            # # # # @show twoμ .* Δϵ̲̲′
            # # # # @show σ̲̲′⁽ᵗʳ⁾
            # # # # @show Si, Hir, rdrsk, κ⁽ᵗʳ⁾, rdrssk, κₛ⁽ᵗʳ⁾
            # # # # # @show d
            # # # # ! update ISV
            # # # ? [20250401T1206] (JMA3): Why are we updating this again?
            # # # * [20250424T0936] (JMA3): Honestly, I have no idea. So let's comment it out for now.
            # # # ! [20250722T1123] (JMA3): It seems as though HEC kept this definition on purpose.
            # # # !                         So let's keep it for now and optimize later.
            # # d       =           (  ( Rs + (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) )  *         ( κₛ ^ (NK-1.0) )          )
            # # # # # @show d
            # # error("Just checking...")
            # if iNewton == 1 # Newton iteration (Backward Euler)
            #     # Nitmax  = 20
            #     # Ntol    = 1.e-06
            #     # Rx0     = X0
            #     # xx      = κ
            #     # for k in range(0, Nitmax)
            #     #     RSRD    = 1.0  +  ( (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) + Rs ) * ψ.Δt * ( xx ^ (NK-1.0) )
            #     #     F1      = (Rx0*κ/RSRD) - xx
            #     #     dF1     = (  -Rx0  *  κ  *  ( (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) + Rs )  *  ψ.Δt  *  ( NK - 1.0 )  *  ( xx ^ (NK-2.0) )  /  ( RSRD ^ 2.0 )  )   -   1.0
            #     #     dxx     = -F1 / dF1
            #     #     xxn     = xx
            #     #     xx     += dxx
            #     #     if abs(dxx / xxn) <= Ntol
            #     #         break
            #     #     end
            #     #     if k >= Nitmax - 1
            #     #         println("Ktr: N-R Conv Issue: k >= Nitmax", k)
            #     #     end
            #     # end
            #     # κ⁽ᵗʳ⁾   = xx
            #     # rdrsk   = 1.0   +   (  ( Rs + (sqrt_twothirds*Rdc*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  ( xx ^ (NK-1.0) )  )
            # end
        #--- trial alpha
            # # @show rs, rd, Δϵ̲̲̇′_mag, ψ.Δt, α̲̲_mag, NK, dzz1
            # error("Just checking...")
            rdrsa   = 1.0   +   (  ( rs + (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  ( α̲̲_mag ^ (NK-1.0) )  *  dzz1  )
            # rdrsa   = 1.0   +   (  ( rs + (                rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  ( α̲̲_mag ^ (NK-1.0) )  *  dzz1  )
            # # @show α̲̲
            # # @show X0, dzz0, rdrsa
            # error("Just checking...")
            α̲̲⁽ᵗʳ⁾ = α̲̲ * (X0*dzz0/rdrsa) # Altr
            # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
            # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
            # # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
            # # rdrsa   = 1.0   -   (  ( rs + (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  ( (sqrt_twothirds*α̲̲_mag) ^ (NK-1.0) )  )
            # # rdrsa   = 1.0   -   (  ( rs + (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  ( (sqrt_threehalves*α̲̲_mag) ^ (NK-1.0) )  *  dzz1  )
            # rdrsa   = 1.0   -   (  ( rs + (sqrt_twothirds* rd*Δϵ̲̲̇′_mag) )  *  ψ.Δt  *  ( (                 α̲̲_mag) ^ (NK-1.0) )  *  dzz1  )
            # α̲̲⁽ᵗʳ⁾ = α̲̲ * (X0*dzz0*rdrsa) # Altr
        #--- Plastic direction tensor N
            ξ̲̲′⁽ᵗʳ⁾      = σ̲̲′⁽ᵗʳ⁾  -  ( (2.0/3.0) .* α̲̲⁽ᵗʳ⁾ ) # Xi
            # ξ̲̲′⁽ᵗʳ⁾      = σ̲̲′⁽ᵗʳ⁾  -  (              α̲̲⁽ᵗʳ⁾ ) # Xi
            # # @show σ̲̲′⁽ᵗʳ⁾
            # # @show α̲̲⁽ᵗʳ⁾
            # error("Just checking...")
            ξ̲̲′⁽ᵗʳ⁾_mag  = norm_symvec(ξ̲̲′⁽ᵗʳ⁾) # Xi_mag
            n̂′          = ξ̲̲′⁽ᵗʳ⁾ ./ ξ̲̲′⁽ᵗʳ⁾_mag # N
            # # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
            # # ξ̲̲′⁽ᵗʳ⁾      = σ̲̲′⁽ᵗʳ⁾  -  ( (2.0/3.0) .* α̲̲⁽ᵗʳ⁾ ) # Xi
            # # ξ̲̲′⁽ᵗʳ⁾_mag  = norm_symvec(ξ̲̲′⁽ᵗʳ⁾) # Xi_mag
            # # n̂′          = ξ̲̲′⁽ᵗʳ⁾ ./ ξ̲̲′⁽ᵗʳ⁾_mag # N
            # # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
            # # ξ̲̲′⁽ᵗʳ⁾      = σ̲̲′⁽ᵗʳ⁾  -  (              α̲̲⁽ᵗʳ⁾ ) # Xi
            # # ξ̲̲′⁽ᵗʳ⁾_mag  = norm_symvec(ξ̲̲′⁽ᵗʳ⁾) # Xi_mag
            # # n̂′          = sqrt_threehalves .* ξ̲̲′⁽ᵗʳ⁾ ./ ξ̲̲′⁽ᵗʳ⁾_mag # N
            # # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
            # ξ̲̲′⁽ᵗʳ⁾      = σ̲̲′⁽ᵗʳ⁾  -  (              α̲̲⁽ᵗʳ⁾ ) # Xi
            # ξ̲̲′⁽ᵗʳ⁾_mag  = norm_symvec(ξ̲̲′⁽ᵗʳ⁾) # Xi_mag
            # n̂′          = ξ̲̲′⁽ᵗʳ⁾ ./ ξ̲̲′⁽ᵗʳ⁾_mag # N
            # # n̂′        ./= norm_symvec(n̂′) # N
    # check plasticity
        ak  = κ⁽ᵗʳ⁾ + Y + β + Yₚ
        ℱ   = (                 ξ̲̲′⁽ᵗʳ⁾_mag) - (sqrt_twothirds*ak*ϕ₁⁽ᵗʳ⁾)
        # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
        # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
        # # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
        # # ℱ   = (sqrt_threehalves*ξ̲̲′⁽ᵗʳ⁾_mag) - (               ak*ϕ₁⁽ᵗʳ⁾)
        # # ℱ   = (                 ξ̲̲′⁽ᵗʳ⁾_mag) - (               ak*ϕ₁⁽ᵗʳ⁾)
        # ℱ   = (                 ξ̲̲′⁽ᵗʳ⁾_mag) - (               ak*ϕ₁⁽ᵗʳ⁾)
    # # # # @show α̲̲_mag, rdrsa
    # # # # @show α̲̲⁽ᵗʳ⁾
    # # # # @show ξ̲̲′⁽ᵗʳ⁾
    # # # # @show ξ̲̲′⁽ᵗʳ⁾_mag
    # # # # @show n̂′
    # # # # @show ak, critra
    # error("Just checking...")
    # Radial-Return
    if ℱ <= 0.0 # elastic solution update
        # @info "Elastic" ℱ σ̲̲′⁽ᵗʳ⁾ α̲̲⁽ᵗʳ⁾ κ⁽ᵗʳ⁾ ϕ η νᵥ ϕ̇ XR XH Xd Xs X
        # deviatoric stress update
        # σ̲̲ = @. σ̲̲⁽ᵗʳ⁾
        σ̲̲′ = σ̲̲′⁽ᵗʳ⁾
        # Cauchy stress update
        σ̲̲ = σ̲̲′ + volumetric(SymmetricTensor{2, 3}(P .* diagm(diag(ones(typeof(σ̲̲)))))) # ! update ISV
        # kinematic hardening & total strain update
        α̲̲ = α̲̲⁽ᵗʳ⁾ # ! update ISV
        # isotropic hardening update
        κ = κ⁽ᵗʳ⁾ # ! update ISV
        # # irradiation hardening update
        # κₛ = κₛ⁽ᵗʳ⁾ # ! update ISV
        # damage update
        ϕ   = ϕ # ! update ISV
        η   = η # ! update ISV
        νᵥ  = νᵥ # ! update ISV
        ϕ̇   = 0.0 # ! update ISV
        # Recrystallization update
        XR  = XR + dXs
        XH  = XH
        Xd  = Xd
        Xs  = Xs + dXs
        X   = X  + dXs
    else # plastic solution (Radial return starts)
        # @info "Plastic" ℱ σ̲̲′⁽ᵗʳ⁾ α̲̲⁽ᵗʳ⁾ κ⁽ᵗʳ⁾ ϕ η νᵥ ϕ̇ XR XH Xd Xs X
        #--- Plastic strain increment solution
        if     iNewton == 0 # analytical solution for DG
            # Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   ak   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
            #     (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   ϕ₁⁽ᵗʳ⁾   *   (  2.0  /  3.0  )   *   (
            #             ( (1.0-X) ^ NK )  *  dzz1  *  ( (h/rdrsa) + (H*Hir/rdrsk) )  )   )    )
            Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   ak   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
                (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   ϕ₁⁽ᵗʳ⁾   *   (  2.0  /  3.0  )   *   (
                        Rx  *  dzz1  *  ( (h/rdrsa) + (H/rdrsk) )  )   )    )
            # # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
            # # Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   (κ⁽ᵗʳ⁾ + Y + (V * asinh(Δϵ̲̲̇′_mag/f)) + Yₚ)   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
            # #     (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   (  2.0  /  3.0  )   *   (
            # #             Rx  *  dzz1  *  ( (h) + (ϕ₁⁽ᵗʳ⁾*H) )  )   )    )
            # # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
            # # Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   (κ⁽ᵗʳ⁾ + Y + (V * asinh(Δϵ̲̲̇′_mag/f)) + Yₚ)   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
            # #     (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   (  2.0  /  3.0  )   *   (
            # #             Rx  *  dzz1  *  ( (h) + (ϕ₁⁽ᵗʳ⁾*H) )  )   )    )
            # # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
            # Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   ak   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
            #     (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   (  2.0  /  3.0  )   *   (
            #             Rx  *  ( (h) + (ϕ₁⁽ᵗʳ⁾*H) )  )   )    )
            # # @show ξ̲̲′⁽ᵗʳ⁾_mag, ak, ϕ₁⁽ᵗʳ⁾, twoμ
            # # @show Rx, dzz1, h, rdrsa, H, rdrsk
            # # @show Δγ
            # error("Just checking...")
            # if t > ψ.Δt
            #     error("Just checking...")
            # end
        elseif iNewton == 1 # Newton-Rapson for DG and Kappa
            # Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   sqrt_twothirds   *   ak   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
            #     (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   ϕ₁⁽ᵗʳ⁾   *   (  2.0  /  3.0  )   *   (
            #             ( (1.0-X) ^ NK )  *  dzz1  *  ( (h/rdrsa) + (H*Hir/rdrsk) )  )   )    )

            # Nitmax  = 20
            # Ntol    = 1e-6
            # xx1     = Δγ
            # xx2     = κ
            # κ₀      = κ
            # Rx      = (1.0-X) ^ NK
            # th      = 1.0 # 1-Backward Euler; 0.5-Midpoint; 0-Forward Euler
            # for k in range(0, Nitmax)
            #     thK0thK     = (  ( 1.0 - th )  *  ( κ₀ ^ (NK-1.0) )  )   +   (  th  *  ( xx2 ^ (NK-1.0) )  )
            #     Rdxx1Rsdt   = 1.0  +  (  ( (sqrt_twothirds*Rdc*xx1) + (Rs*ψ.Δt) )  *  thK0thK  )
            #     F₁      = ξ̲̲′⁽ᵗʳ⁾_mag   -   (  twoμ * xx1  )   -   (
            #             sqrt_twothirds  *  ( κ₀ + (Rx*H*Hir*xx1) )  /  Rdxx1Rsdt  )   -   (
            #             sqrt_twothirds  *  ( Be + Y + Yₚ )  )
            #     ∂F₁╱∂x₁ = -twoμ     -     (    sqrt_twothirds    *    (
            #         (  Rx  *  H  *  Hir  /  Rdxx1Rsdt  )   -   (
            #             ( κ₀ + (Rx*H*Hir*xx1) )  *  sqrt_twothirds  *  Rdc  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )   )    )
            #     ∂F₁╱∂x₂ = -1.0 * sqrt_twothirds
            #     F₂      = (  ( κ₀ + (Rx*H*Hir*xx1) )  /  Rdxx1Rsdt  )   -   xx2
            #     ∂F₂╱∂x₁ = (  Rx  *  H  *  Hir  /  Rdxx1Rsdt  )   -   (
            #         ( κ₀ + (Rx*H*Hir*xx1) )  *  sqrt_twothirds  *  Rdc  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )
            #     ∂F₂╱∂x₂ = -( κ₀ + (Rx*H*Hir*xx1) )  *  (
            #         (sqrt_twothirds*Rdc*xx1) + (Rs*ψ.Δt) )  *  th  *  ( NK - 1.0 )  *  ( xx2 ^ (NK-2.0) )
            #     ∂F₂╱∂x₂ = ( ∂F₂╱∂x₂ / (Rdxx1Rsdt^2.0) )  -  1.0

            #     a₁₁     = ∂F₁╱∂x₁
            #     a₁₂     = ∂F₁╱∂x₂
            #     a₂₁     = ∂F₂╱∂x₁
            #     a₂₂     = ∂F₂╱∂x₂
            #     dxx2    = ( (-F₂*a₁₁/a₂₁) + F₁ )  /  ( (a₂₂ * (a₁₁/a₂₁)) - a₁₂ )
            #     dxx1    = ( (-a₁₂*dxx2)   - F₁ )  /  a₁₁

            #     xx1n    = xx1
            #     xx2n    = xx2

            #     xx1    += dxx1
            #     xx2    += dxx2

            #     if abs(dxx1/xx1n) <= Ntol && abs(dxx2/xx2n) <= Ntol
            #         break
            #     end
            #     if k >= Nitmax - 1
            #         println("Gamma-K: N-R Conv. Issue: k >= Nitmax")
            #     end
            # end
            # Δγ  = xx1
            # κ   = xx2 # ! update ISV



            Δγ = (    ξ̲̲′⁽ᵗʳ⁾_mag    -    (   ak   *   ϕ₁⁽ᵗʳ⁾   )    )     /     (
                (   ϕ₁⁽ᵗʳ⁾   *   twoμ   )    +    (   ϕ₁⁽ᵗʳ⁾   *   (  2.0  /  3.0  )   *   (
                        Rx  *  dzz1  *  ( (h/rdrsa) + (H/rdrsk) )  )   )    )

            Nitmax  = 100
            Ntol    = 1e-4
            xx1     = Δγ
            xx2     = κ
            κ₀      = κ
            # Rx      = (1.0-X) ^ NK
            th      = 1.0 # 1-Backward Euler; 0.5-Midpoint; 0-Forward Euler
            for k in range(0, Nitmax)
                thK0thK     = (  ( 1.0 - th )  *  ( κ₀ ^ (NK-1.0) )  )   +   (  th  *  ( xx2 ^ (NK-1.0) )  )
                Rdxx1Rsdt   = 1.0  +  (  ( (sqrt_twothirds* Rd*xx1) + (Rs*ψ.Δt) )  *  thK0thK  )
                F₁      = ξ̲̲′⁽ᵗʳ⁾_mag   -   (  twoμ * ϕ₁⁽ᵗʳ⁾ * xx1  )   -   (
                        sqrt_twothirds  *  ( κ₀ + (Rx*H*ϕ₁⁽ᵗʳ⁾*xx1) )  /  Rdxx1Rsdt  )   -   (
                        sqrt_twothirds  *  ( β + Y + Yₚ )  )
                ∂F₁╱∂x₁ = -twoμ     *     ϕ₁⁽ᵗʳ⁾     -     (    sqrt_twothirds    *    (
                    (  Rx  *  H  *  ϕ₁⁽ᵗʳ⁾  /  Rdxx1Rsdt  )   -   (
                        ( κ₀ + (Rx*H*ϕ₁⁽ᵗʳ⁾*xx1) )  *  sqrt_twothirds  *  Rd  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )   )    )
                ∂F₁╱∂x₂ = -1.0 * sqrt_twothirds
                F₂      = (  ( κ₀ + (Rx*H*ϕ₁⁽ᵗʳ⁾*xx1) )  /  Rdxx1Rsdt  )   -   xx2
                ∂F₂╱∂x₁ = (  Rx  *  H  *  ϕ₁⁽ᵗʳ⁾  /  Rdxx1Rsdt  )   -   (
                    ( κ₀ + (Rx*H*ϕ₁⁽ᵗʳ⁾*xx1) )  *  sqrt_twothirds  *  Rd  *  thK0thK  /  ( Rdxx1Rsdt ^ 2.0 )  )
                ∂F₂╱∂x₂ = -( κ₀ + (Rx*H*ϕ₁⁽ᵗʳ⁾*xx1) )  *  (
                    (sqrt_twothirds*Rd*xx1) + (Rs*ψ.Δt) )  *  th  *  ( NK - 1.0 )  *  ( xx2 ^ (NK-2.0) )
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
                    println("Gamma-K: N-R Conv. Issue: k >= Nitmax", abs(dxx2/xx2n))
                end
            end
            Δγ  = xx1
            κ   = xx2 # ! update ISV
        end
        #--- stress solution
            # deviatoric stress update
            # # @show σ̲̲′⁽ᵗʳ⁾
            # # @show ϕ₁⁽ᵗʳ⁾, twoμ, Δγ, ϕ₁⁽ᵗʳ⁾*twoμ*Δγ
            # # @show n̂′
            # # @show σ̲̲′⁽ᵗʳ⁾  -  ( (ϕ₁⁽ᵗʳ⁾*twoμ*Δγ) .* n̂′ )
            σ̲̲′ = σ̲̲′⁽ᵗʳ⁾  -  ( (ϕ₁⁽ᵗʳ⁾*twoμ*Δγ) .* n̂′ )
            # Cauchy stress update
            # σ̲̲ = σ̲̲′ + volumetric(P) # ! update state variable
            # @show σ̲̲
            # @show hydrostatic(σ̲̲)
            # @show ϕ₁⁽ᵗʳ⁾, ψ.Δt, K
            # @show ϵ̲̲
            # @show I₁(ϵ̲̲)
            # @show hydrostatic(σ̲̲)*ϕ₁⁽ᵗʳ⁾
            # @show ϕ₁⁽ᵗʳ⁾*ψ.Δt*K*I₁(ϵ̲̲)
            # @show (hydrostatic(σ̲̲)*ϕ₁⁽ᵗʳ⁾) + (ϕ₁⁽ᵗʳ⁾*ψ.Δt*K*I₁(ϵ̲̲))
            # @show ((hydrostatic(σ̲̲)*ϕ₁⁽ᵗʳ⁾) + (ϕ₁⁽ᵗʳ⁾*ψ.Δt*K*I₁(ϵ̲̲))) .* diagm(diag(ones(typeof(σ̲̲))))
            # @show SymmetricTensor{2, 3}(((hydrostatic(σ̲̲)*ϕ₁⁽ᵗʳ⁾) + (ϕ₁⁽ᵗʳ⁾*ψ.Δt*K*I₁(ϵ̲̲))) .* diagm(diag(ones(typeof(σ̲̲)))))
            # @show σ̲̲′
            σ̲̲ = σ̲̲′ + SymmetricTensor{2, 3}(((hydrostatic(σ̲̲)*ϕ₁⁽ᵗʳ⁾) + (ϕ₁⁽ᵗʳ⁾*ψ.Δt*K*I₁(ϵ̲̲))) .* diagm(diag(ones(typeof(σ̲̲))))) # ! update state variable
            # # # @show σ̲̲′⁽ᵗʳ⁾
            # # # @show ϕ₁⁽ᵗʳ⁾, twoμ, Δγ
            # # # @show n̂′
            # # # # @show σ̲̲′⁽ᵗʳ⁾  -  ( (ϕ₁⁽ᵗʳ⁾*(89066.17)*(0.01089263427151211)) .* n̂′ )
            # # # @show σ̲̲′, vM
            # # # @show σ̲̲
            # # error("Just checking...")
            # if t > ψ.Δt
            #     error("Just checking...")
            # end
        #--- total plastic strain
            ϵ̲̲⁽ᵖ⁾ += ( (sqrt_twothirds*Δγ) .* n̂′ ) # ! update ISV
            # ϵ̲̲⁽ᵖ⁾ += ( (               Δγ) .* (n̂′) ) # ! update ISV
            # # # @show norm_symvec(ϵ̲̲⁽ᵖ⁾)
        #--- alpha solution
            # α̲̲ = α̲̲⁽ᵗʳ⁾    +    ( # ! update ISV
            #         (  ( (1.0-X) ^ NK )  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  h  *  Δγ  )   .*   n̂′   ./   rdrsa   )
            α̲̲ = α̲̲⁽ᵗʳ⁾    +    ( # ! update ISV
                    (  Rx  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  h                     *  Δγ  )   .*   n̂′   ./   rdrsa   )
            # α̲̲ = α̲̲⁽ᵗʳ⁾    +    ( # ! update ISV
            #         (  Rx  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  h  *  sqrt_twothirds  *  Δγ  )   .*   n̂′   ./   rdrsa   )
            # α̲̲ = α̲̲⁽ᵗʳ⁾    +    ( # ! update ISV
            #         (  Rx  *  dzz1  *  h  *  Δγ  )   .*   n̂′   )
            # α̲̲ = α̲̲⁽ᵗʳ⁾    +    ( # ! update ISV
            #         (  Rx  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  h  *  Δγ  )   .*   n̂′   )
        #--- kappa solution # ! update ISV
            # κ = (   iNewton   !=   0   )    ?    (   xx2   )    :    (#=[=#   κ⁽ᵗʳ⁾   +   (
            #         ( (1.0-X) ^ NK )  *  sqrt_twothirds  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  H  *  Hir  *  Δγ  /  rdrsk  )   #=]=#)
            κ = (   iNewton   !=   0   )    ?    (   xx2   )    :    (#=[=#   κ⁽ᵗʳ⁾   +   (
                    Rx  *  sqrt_twothirds  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  H  *  Δγ  /  rdrsk  )   #=]=#)
            # κ = (   iNewton   !=   0   )    ?    (   xx2   )    :    (#=[=#   κ⁽ᵗʳ⁾   +   (
            #         Rx                     *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  H  *  Δγ  /  rdrsk  )   #=]=#)
            # κ = (   iNewton   !=   0   )    ?    (   xx2   )    :    (#=[=#   κ⁽ᵗʳ⁾   +   (
            #         Rx  *  dzz1  *  sqrt_twothirds  *  H  *  Δγ  )   #=]=#)
            # κ = (   iNewton   !=   0   )    ?    (   xx2   )    :    (#=[=#   κ⁽ᵗʳ⁾   +   (
            #         Rx  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  sqrt_twothirds  *  H  *  Δγ  )   #=]=#)
        # #--- irradiation hardening solution in isotropic hardening
            # Sir = Si^2.0
            # κₛ  = κₛ⁽ᵗʳ⁾   +   ( # ! update ISV
            #     ( (1.0-X) ^ NK )  *  sqrt_twothirds  *  dzz1  *  ϕ₁⁽ᵗʳ⁾  *  H  *  Sir  *  Δγ  /  rdrssk  )
        # various print statements for debugging
            # # # # @show Δγ
            # # # # @show σ̲̲′
            # # # # @show σ̲̲, vM
            # # # # @show ϵ̲̲′
            # # # # @show ϵ̲̲⁽ᵖ⁾
            # # # # @show ϵ̲̲⁽ᴴ⁾
            # # # # @show α̲̲
            # # # # @show κ
            # error("Just checking...")
        #--- damage
            di1 = I₁(σ̲̲)
            dj2 = I₂(σ̲̲′)
            dj3 = I₃(σ̲̲′)
            # # # # @show σ̲̲′, di1, dj2, dj3
            JJ1 = (dj3^2.0) / (dj2^3.0)
            JJ2 = (dj3    ) / (dj2^1.5)
            JJ3 = (di1    ) / (dj2^0.5)
            # # [20250401T1042] (JMA3): this comment (v) is from HEC's original code
            # # here I controlled stress triaxiality to 1 (tension (Horstemeyer et al., 2000))
            # JJ3 = 1.0
            # # # # @show σ̲̲′
            # # # # @show di1, dj2, dj3, JJ1, JJ2, JJ3
        ##--- nucleation (RK4 integration)
            ddff= ( ψ.𝒹 ^ 0.5 )  /  ( ψ.𝒻 ^ (1.0/3.0) )
            # # Δη₀ = Δϵ̲̲̇′_mag   *   ddff   /   ψ.Kic   *   (  a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (
            # #     c * damirr * abs(JJ3) )  )   *   exp(  Tnuc  /  ψ.θ  )
            # #     #+ pcc*(1.+sinh(kp1*Si))*abs(JJ3))*exp(pTnuc/ψ.θ)
            # η̇₀  = Δϵ̲̲̇′_mag   *   ddff   /   ψ.Kic   *   (  a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (
            #     c * abs(JJ3) )  )   *   exp(  Tnuc  /  ψ.θ  )
            #     #+ pcc*(1.+sinh(kp1*Si))*abs(JJ3))*exp(pTnuc/ψ.θ)
            # k₁  = η̇₀  *    η
            # k₂  = η̇₀  *  ( η + (0.5k₁*ψ.Δt) )
            # k₃  = η̇₀  *  ( η + (0.5k₂*ψ.Δt) )
            # k₄  = η̇₀  *  ( η + (   k₃*ψ.Δt) )
            # # # # # @show η
            # η  += 6.0  \  ψ.Δt  *  ( k₁ + 2.0(k₂+k₃) + k₄ ) # ! update ISV
            # # η̇  = η   *   Δϵ̲̲̇′_mag   *   ddff   /   ψ.Kic   *   (  a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (
            # #     c * damirr * abs(JJ3) )  )   *   exp(  Tnuc  /  ψ.θ  )
            # η̇   = η   *   Δϵ̲̲̇′_mag   *   ddff   /   ψ.Kic   *   (  a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (
            #     c * abs(JJ3) )  )   *   exp(  Tnuc  /  ψ.θ  )
            # # # # # # @show k₁, k₂, k₃, k₄
            # # # # # @show ddff, Δη₀, η, Δη

            ### Implementation (Horstemeyer et al., 2000)
            # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
            # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
            η₀  = η
            # η   = Cnuc * exp(   ( ϵ̲̲′_mag + (Δϵ̲̲̇′_mag*ψ.Δt) )   *   ddff   /   ψ.Kic   *   (
            #     a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (c * abs(JJ3) )  )   )
            # η   = Cnuc * exp(   ( ϵ̲̲′_mag + (Δϵ̲̲̇′_mag*ψ.Δt) )   *   ddff   /   ψ.Kic   *   (
            #     a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (c * abs(JJ3) )  )   )    *    exp(   Tnuc   /   ψ.θ   )
            η   = Cnuc * exp(   ϵ̲̲′_mag   *   ddff   /   ψ.Kic   *   (
                a  *  ( (4.0/27.0) - JJ1 )  +  ( b * JJ2 )  +  (c * abs(JJ3) )  )   )    *    exp(   Tnuc   /   ψ.θ   )
            η̇   = (η-ψ.η₀) / ψ.Δt
            # η̇   = (η-η₀) / ψ.Δt

            #nuc0 = PE[i]*ddff/pKic*(paa*(4./27.-JJ1) + pbb*(JJ2) \
            #     + pcc*(1.+sinh(pccsi*Si))*abs(JJ3))*exp(pTnuc/ψ.θ)
            #Nuc[i] = Cnuc*exp(nuc0)
        ##--- growth
            # ### Implementation (Euler method)
            # #dvod = 4./3.*((sqrt(3.)/2.*prr0*ddd/(1.-pnn) \
            # #     * sinh(sqrt(3.)*(1.-pnn)*sqrt(2.)/3.*JJ3)) \
            # #     * exp(pTgrw*temp))**3
            # ν̇ᵥ  = (  ( √(3.0) / 2.0 )  *  ψ.R₀  *  ϵ̲̲′_mag  /  ( 1.0 - nn )
            #         * sinh( √(3.0) * (1.0-nn) * (√(2.0)/3.0) * JJ3 )
            #     )   *   exp(  Tgrw  *  ψ.θ  )   *   νᵥ
            # νᵥ += ν̇ᵥ*ψ.Δt

            ### Implementation (Horstemeyer et al., 2000)
            # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
            # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
            # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
            νᵥ₀ = νᵥ
            # ! update ISV
            # νᵥ = (    4.0    /    3.0    )     *     (#={=#    (   ψ.R₀   *   exp(#=[=#
            #         ( ϵ̲̲′_mag + (Δϵ̲̲̇′_mag*ψ.Δt) )  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-nn) )  *  sinh(
            #             sqrt(3.0) * (1.0-nn) * sqrt(2.0) / 3.0 * JJ3 )
            #     #=]=#)   )    ^    3.0    #=}=#)
            # νᵥ = ψ.R₀   *   (#=[=#
            #         Δϵ̲̲̇′_mag  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-nn) )  *  sinh(
            #             sqrt(3.0) * (1.0-nn) * sqrt(2.0) / 3.0 * JJ3 )
            #     #=]=#)
            # νᵥ = (    4.0    /    3.0    )     *     (#={=#    (   ψ.R₀   *   exp(#=[=#
            #         ( ϵ̲̲′_mag + (Δϵ̲̲̇′_mag*ψ.Δt) )  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-nn) )  *  sinh(
            #             sqrt(3.0) * (1.0-nn) * sqrt(2.0) / 3.0 * JJ3 )  *  exp( Tgrw * ψ.θ )
            #     #=]=#)   )    ^    3.0    #=}=#)
            νᵥ = (    4.0    /    3.0    )     *     (#={=#    (   ψ.R₀   *   begin
                    x = exp(#=[=#
                        ϵ̲̲′_mag  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-nn) )  *  sinh(
                            sqrt(3.0) * (1.0-nn) * sqrt(2.0) / 3.0 * JJ3 )  *  exp( Tgrw * ψ.θ )
                    #=]=#)
                    x = isinf(x) ? 0.0 : x
                end)    ^    3.0    #=}=#)
            ν̇ᵥ = (νᵥ-νᵥ₀) / ψ.Δt
            # @show νᵥ₀, ψ.R₀, ϵ̲̲′_mag, nn, JJ3, Tgrw, ψ.θ, νᵥ, ν̇ᵥ
            # @show ϵ̲̲′_mag  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-nn) )
            # @show sinh(sqrt(3.0) * (1.0-nn) * sqrt(2.0) / 3.0 * JJ3 )
            # @show exp( Tgrw * ψ.θ )
            # @show (   exp(#=[=# ϵ̲̲′_mag  *  sqrt( 3.0 )  /  ( 2.0 * (1.0-nn) )  *  sinh( sqrt(3.0) * (1.0-nn) * sqrt(2.0) / 3.0 * JJ3 )  *  exp( Tgrw * ψ.θ ) #=]=#)   )

            # # # # @show νᵥ₀, ϵ̲̲′_mag, νᵥ, Δνᵥ

            #vod0 = PE[i]*sqrt(3.)/(2.*(1.-pnn)) \
            #     * sinh(sqrt(3.)*(1.-pnn)*sqrt(2.)/3.*JJ3)*exp(pTgrw*ψ.θ)
            #Vod[i] = 4./3.*(prr0*exp(vod0))^3
        ##--- coalesence
        C = 1.0 # ! update ISV
        ##--- damage rate
        ϕ̇ = (η̇*νᵥ) + (η*ν̇ᵥ) # ! update ISV
        # # # # @show ϕ̇
        ##--- total damage at current step
        # ? [2025T1048] (JMA3): why the blazes does this phi have 3 re-assignments?
        ϕ₀    = ϕ
        # # # # @show ϕ₀
        # # Horstemeyer, et. al. (2000): https://www.sciencedirect.com/science/article/pii/S016784429900049X?ref=pdf_download&fr=RR-2&rr=9674f0426d71c595
        # # modified in Cho, et. al. (2017): https://onlinelibrary.wiley.com/doi/epdf/10.1002/9781119018377.ch7?saml_referrer=
        # # modified (again) in Cho, et. al. (2019): https://www.sciencedirect.com/science/article/pii/S0749641918303139?casa_token=tQbSk0wbfLwAAAAA:vQJyOp3-HPScV3EmVpZOT3Hpx6cCBa_Gwft4WzdFHHLRqSpD1s66BdkpqM8BIl4AC-Qn1bUDZg
        # # ϕ   = C*η*νᵥ
        # # ϕ   = C*η*νᵥ * ((ψ.d₀/d)^ψ.z) * exp(Tgrw*ψ.θ)
        # ϕ   = C*η*νᵥ
        # # # # @show ϕ
        ϕ  += ϕ̇*ψ.Δt
        # # # # @show ϕ
        # ϕ     = max( min(ϕ,0.99999) , 0.0000000001 ) # ! update ISV
        ϕ   = max( min(ϕ,0.999999999) , 0.0000000001 ) # ! update ISV
        # # # # @show ϕ

        ϕ̇ = (ϕ-ϕ₀) / ψ.Δt # ! update ISV

        # # # # @show ϕ₀, η, νᵥ, C, ϕ, ϕ̇, ϵ̲̲⁽ᴴ⁾
        # error("Just checking...")
    end
    # # @show t
    # # @show ϕ₁⁽ᵗʳ⁾
    # # @show twoμ
    # # @show Δγ
    # # @show n̂′
    # # @show σ̲̲
    # # @show deviatoric(σ̲̲)
    # # @show vM
    # # @show ϵ̲̲
    # # @show deviatoric(ϵ̲̲)
    # # @show ϵ̲̲⁽ᵖ⁾
    # # # @show # # @show ϵ̲̲⁽ᵖ⁾[1] + (sqrt_twothirds*Δγ)
    # # @show # # @show norm_symvec(ϵ̲̲⁽ᵖ⁾)
    # # @show α̲̲
    # # @show κ
    # # @show ϕ
    # # @show η
    # # @show νᵥ
    # # @show ϕ̇
    # # @show X
    # # @show d
    # # error("Just checking...")
    # if t > 50Δt
    #     error("Just checking...")
    # end
    return MaterialState(SymmetricTensor{2, 3}(σ̲̲), ϵ̲̲, SymmetricTensor{2, 3}(ϵ̲̲⁽ᵖ⁾), SymmetricTensor{2, 3}(α̲̲), κ, ϕ, η, νᵥ, ϕ̇, X, XR, XH, Xs, Xd, d)
end

function ContinuumMechanicsBase.predict(
            ψ   ::Cho2019UnifiedStaticDynamicTensor{T, S},
            test::BammannChiesaJohnsonPlasticity.AbstractBCJMetalTest{T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat, S<:SymmetricTensor{4, 3, T}}
    M = ψ.N + 1
    # numerical constants
    sqrt_threehalves = √(3.0/2.0)
    # irradiation before damage
    Tirr = ψ.P
    M0, Si, damirr = 0.0, 0.0, 1.0
    if Tirr != 0.0
        kr = kr1 * exp(krt/Tirr)
        Si = (kr*flu) ^ (1.0/kr2)
        M0 = Kr3 * Si
        damirr = exp(  ( kp1 * exp(kpt/Tirr) * flu )  ^  ( 1.0 / kp2 )  )
    end


    ###########################################################
    # define RVE
    L = 1.0
    grid = generate_grid(Hexahedron, (1, 1, 1), Vec{3}([0., 0., 0.]), Vec{3}([1., 1., 1.] .* L))
    addfacetset!(grid, "back_yz", x-> x[1] ≈ 0.)
    addfacetset!(grid, "back_xz", x-> x[2] ≈ 0.)
    addfacetset!(grid, "back_xy", x-> x[3] ≈ 0.)
    addfacetset!(grid, "front_xz", x-> x[2] ≈ L)
    # addcellset!(grid, "front_xz", x-> x[2] ≈ L)
    addnodeset!(grid, "rve", x->x[1] >= 0.)
    addnodeset!(grid, "front_xz", x->x[2] ≈ L)
    gridnodeset = getnodeset(grid, "rve")
    refshape = RefHexahedron
    gridncells = getncells(grid)
    gridnnodes = getnnodes(grid)
    interpolation = Lagrange{refshape, 1}()^3
    ## degrees of freedom
    # dh = create_dofhandler(grid, interpolation) # JuaFEM helper function
    dh = DofHandler(grid)
    add!(dh, :u, interpolation) # add a displacement field with 3 components
    close!(dh)

    vectorstrafe_i(i) = vectorstrafe(grid, dh, i)

    ## constraints
    # dbcs = create_bc(dh, grid) # create Dirichlet boundary-conditions
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(grid, "back_yz"), (x, t) -> [0.0], [1]))
    add!(ch, Dirichlet(:u, getfacetset(grid, "back_xz"), (x, t) -> [0.0], [2]))
    add!(ch, Dirichlet(:u, getfacetset(grid, "back_xy"), (x, t) -> [0.0], [3]))
    close!(ch)

    # define boundary conditions and time steps
    # displacement = do_verify ? collect(range(first(expdata_strain), last(expdata_strain); length=n_timesteps)) : nothing
    displacement = collect(range(0.0, ψ.ϵₙ; length=M))
    # V = -1
    # V = do_verify ? expdata_ϵ̇ : -1.3477876909913353 # initial velocity onto Right plate in local reference frame
    V = ψ.ϵ̇_eff # initial velocity onto Right plate in local reference frame
    # time_0 = do_verify ? (first(expdata_strain) / V) : 0.0
    # time_n = do_verify ? (last(expdata_strain) / V) : 0.1
    time_0 = 0.0
    time_n = ψ.ϵₙ / V
    time_domain = collect(range(time_0, time_n; length=M))
    VV = zeros(M)
    ΔX = zeros(M)

    # pre-allocate solution vectors, etc.
    # ! ndofs() ≡ dim * getnnodes(grid)
    n_dofs  = ndofs(dh)  # total number of dofs
    u       = zeros(n_dofs)  # solution vector
    Δu      = zeros(n_dofs)  # displacement correction
    u₀      = zeros(n_dofs)
    r       = zeros(n_dofs)   # residual
    K       = allocate_matrix(dh, ch); # tangent stiffness matrix

    # create material states. One array for each cell, where each element is an array of material-
    # states - one for each integration point
    ## quadratures
    # cellvalues_u, facetvalues_u = create_values(interpolation)
    # # cellvalues_u, cellvalues_eps, facetvalues = create_values(interpolation)
    # setup quadrature rules
    qr      = QuadratureRule{refshape}(2)
    facet_qr = FacetQuadratureRule{refshape}(refshape == RefHexahedron ? 2 : 3)
    # cell and facetvalues for u
    cellvalues = CellValues(qr, interpolation)
    facetvalues = FacetValues(facet_qr, interpolation)
    nu = getnbasefunctions(cellvalues)
    cell_basefuncs = getnbasefunctions(cellvalues)
    facet_basefuncs = getnbasefunctions(facetvalues)

    ## material states
    nqp_cell = getnquadpoints(cellvalues)
    nqp_facet = getnquadpoints(facetvalues)
    # TODO [20250214T0935] (JMA3) - figure out how to restart from previous calculation in case of interruptions
    # states_material = [MaterialState() for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_material_old = [MaterialState() for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_material = [JohnsonCookState() for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_material_old = [JohnsonCookState() for _ in 1:nqp_cell, _ in 1:gridncells]
    states_material = [[MaterialState(ψ) for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:M]
    # states_strain = [[zero(SymmetricTensor{2, 3}) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells]
    states_strain = [[[zero(SymmetricTensor{2, 3}) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:M]
    # states_strain_old = [[zero(SymmetricTensor{2, 3}) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_strain_rate = [[zero(SymmetricTensor{2, 3}) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells]
    states_strain_rate = [[[zero(SymmetricTensor{2, 3}) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:M]
    # states_energy_total = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_total_old = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_elastic = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_elastic_old = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_plastic = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_plastic_old = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    states_energy_total = [[0. for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:M]
    states_energy_elastic = [[0. for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:M]
    states_energy_plastic = [[0. for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:M]
    # states_traction = [zeros((3, 3)) for _ in 1:nqp_cell, _ in 1:gridncells]
    states_traction = [[zeros((3, 3)) for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:M]
    # states_traction_old = [zeros((3, 3)) for _ in 1:nqp_cell, _ in 1:gridncells]

    # create export vectors results file
    u_max           = zeros(M)
    strain_max      = zeros(M)
    strainrate_max  = zeros(M)
    stress_max      = zeros(M)
    α_max           = zeros(M) # hcat(α⃗...)
    κ_max           = zeros(M) # hcat(κ_vec...)
    ϕ_max           = zeros(M) # hcat(ϕ⃗...)
    X_max           = zeros(M) # hcat(X⃗...)
    d_max           = zeros(M) # hcat(d⃗...)
    pvd             = paraview_collection("rve")

    gridnodes = deepcopy(getnodes(grid))
    gridnodes_view = @view gridnodes[[1:gridnnodes...]]
    function update_rhs(grid, n)
        dx = get_node_coordinate(grid, n) - get_node_coordinate(gridnodes[node_sphere_o])
        Δx = (R * (dx / norm(dx))) - dx
        return Δx
    end
    function update_rhs(grid, n, v)
        px, py, pz = get_node_coordinate(gridnodes[n])
        cx, cy, cz = get_node_coordinate(gridnodes[node_sphere_o])
        vx, vy, vz = v
        # Calculate A, B, C for the quadratic equation
        A = vx^2 + vy^2 + vz^2
        B = 2 * ((px - cx) * vx + (py - cy) * vy + (pz - cz) * vz)
        C = (px - cx)^2 + (py - cy)^2 + (pz - cz)^2 - R^2

        # Discriminant of the quadratic equation
        discriminant = B^2 - 4A * C
        if discriminant < 0
            error("No real intersection: direction vector may be incorrect.")
        end

        # Solve for t using the quadratic formula
        t1 = (-B + sqrt(discriminant)) / (2A)
        t2 = (-B - sqrt(discriminant)) / (2A)

        # Choose the smallest positive t
        t_min = min(filter(x -> x > 0, [t1, t2])...)

        # Calculate the vector from P to the surface (Q)
        vector_to_surface = t_min * [vx, vy, vz]

        return vector_to_surface
    end
    function update_rhs_view(grid, n)
        return get_node_coordinate(gridnodes[node_sphere_z]) - get_node_coordinate(gridnodes[n])
    end
    ###########################################################



    # # begin prediction
    # σ⃗ = []; push!(σ⃗, σ̲̲)
    # ϵ⃗ = []; push!(ϵ⃗, ϵ̲̲)
    # α⃗ = []; push!(α⃗, 0.0)
    # κ_vec = []; push!(κ_vec, 0.0)
    # ϕ⃗ = []; push!(ϕ⃗, ϕ)
    # X⃗ = []; push!(X⃗, X)
    # d⃗ = []; push!(d⃗, d)
    # ϵ̲̲       = zero(SymmetricTensor{2, 3}) # zeros(ψ.Δϵ̲̲)
    t       = 0.0

    # march through time
    NEWTON_TOL = 1e9 # 1 # 1 N
    NEWTON_M = 10
    # println("\nStarting Netwon iterations:")
    # Δt = do_verify ? (time_domain[2] - time_domain[1]) : (time_n / (n_timesteps - 1))
    for i ∈ range(2, M)
    # for i ∈ range(2, ψ.N)
    # for i ∈ range(2, 3)
        t += ψ.Δt
        # ϵ̲̲ += ψ.Δϵ̲̲
        # if i > 3
        #     error("Just checking...")
        # end
        ###########################################################
        # @printf("\n Time step %d (t = %.6f s) @ %.4f m/s:\n", timestep - 1, t, V) # -1 to match ParaView
        VV[i] = V
        ΔV = 0.
        # p = 0.
        # state_material = states_material
        # state_material_old = states_material_old
        # # println(first(state_material).σ == first(state_material_old).σ)
        state_material_old = (i > 1 ? states_material[i - 1] : states_material[1])
        state_material = (i > 1 ? states_material[i - 1] : states_material[1])
        # state_strain = states_strain
        # state_strain_old = states_strain_old
        state_strain_old = (i > 1 ? states_strain[i - 1] : states_strain[1])
        state_strain = (i > 1 ? states_strain[i - 1] : states_strain[1])
        state_strain_itr = deepcopy(state_strain)
        # state_energy_total = states_energy_total
        # state_energy_elastic = states_energy_elastic
        # state_energy_plastic = states_energy_plastic
        # state_energy_total_old = states_energy_total_old
        # state_energy_total_itr = deepcopy(state_energy_total_old)
        state_energy_total = states_energy_total[i]
        state_energy_elastic = states_energy_elastic[i]
        state_energy_plastic = states_energy_plastic[i]
        state_energy_total_old = (i > 1 ? states_energy_total[i - 1] : state_energy_total)
        state_energy_total_itr = deepcopy(state_energy_total_old)
        # state_traction = states_traction
        # state_traction_old = states_traction_old
        state_traction = states_traction[i]
        state_traction_old = (i > 1 ? states_traction[i - 1] : states_traction[1])
        # # (v)    only needed if ch[f(x, t)]     (v)
        # update!(ch, t) # evaluates the D-bndc at time t
        # # (^)                                   (^)
        apply!(u, ch)  # set the prescribed values in the solution vector
        if i > 1 # t > 0
            # ΔX[i] = ψ.Δϵ̲̲
            ΔX[i] = displacement[i]
            dbcs_t = ConstraintHandler(dh)
            add!(dbcs_t, Dirichlet(:u, getfacetset(grid, "front_xz"), (x, t) -> [ΔX[i]], [2]))
            close!(dbcs_t)
            apply!(u, dbcs_t)  # set the prescribed values in the solution vector
        end
        # @printf("\n Time step %d (t = %.6f s) @ %.6f mm:\n", i - 1, t, ΔX[i]) # -1 to match ParaView
        # @show t, ΔX[i]
        println_i_rve(n) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m]\n",
            n, get_node_coordinate(gridnodes[n])[2], u[vectorstrafe_i(n)][2], ΔX[n])
        println_i_rve(n, trac, et, ee, ep, p, ΔV) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m] | T:%.3e [Pa] | Eᵗₘₐₓ:%.3e [J] | Eᵉₘₐₓ:%.3e [J] | Eᵖₘₐₓ:%.3e [J] | %%_diff(E):%.3f [%%] | ΔV:%.6f [m/s]\n",
            n, get_node_coordinate(gridnodes[n])[2], u[vectorstrafe_i(n)][2], ΔX[n], maximum(abs, map(norm, trac)), maximum(abs, et), maximum(abs, ee), maximum(abs, ep), p, ΔV)

        # newton-raphson loop
        newton_itr = 1
        while newton_itr <= NEWTON_M
            # Tangent and residual contribution from the cells (volume integral)
            # doassemble!(K, r, cellvalues, dh, material, u, states, states_old);
            assembler = start_assemble(K, r)
            # assembler = start_assemble(K, r; fillzero=false)
            if t > 0.
                re = zeros(nu)     # element residual vector
                ke = zeros(nu, nu) # element tangent matrix
                for cell in CellIterator(dh)
                    fill!(ke, 0)
                    fill!(re, 0)
                    eldofs = celldofs(cell)
                    ue = u[eldofs]
                    # @show ue
                    # re = @view r[eldofs]
                    # ke = @view K[eldofs, eldofs]
                    state = @view state_material[:, cell.cellid]
                    state_old = @view state_material_old[:, cell.cellid]
                    strain = @view state_strain[:, cell.cellid]
                    strain_itr = @view state_strain_itr[:, cell.cellid]
                    strain_old = @view state_strain_old[:, cell.cellid]
                    material = ψ
                    stateenergy_total = @view state_energy_total[:, cell.cellid]
                    stateenergy_elastic = @view state_energy_elastic[:, cell.cellid]
                    stateenergy_plastic = @view state_energy_plastic[:, cell.cellid]
                    statetraction = @view state_traction[:, cell.cellid]
                    fill!(stateenergy_total,    0.0)
                    fill!(stateenergy_elastic,  0.0)
                    fill!(stateenergy_plastic,  0.0)
                    # fill!(statetraction,        0.0)
                    Ferrite.reinit!(cellvalues, cell)
                    for q_point in 1:nqp_cell
                        # for each integration point, compute stress and material stiffness
                        # ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
                        # σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state_old[q_point])
                        Δϵ = function_symmetric_gradient(cellvalues, q_point, ue) - state[q_point].ϵ̲̲ # Total strain
                        # @show state[q_point].ϵ̲̲
                        # @show Δϵ
                        # @show state[q_point].ϵ̲̲ + Δϵ
                        σ_prev, ϵ_prev, ϵᵖ_prev = state[q_point].σ̲̲, state[q_point].ϵ̲̲, state[q_point].ϵ̲̲⁽ᵖ⁾
                        # if q_point == 1
                        #     # println(size(function_symmetric_gradient(cellvalues, q_point, ue)))
                        #     println(size(shape_gradient(cellvalues, q_point, 1)))
                        # end
                        # Δϵ, σ, D, state[q_point] = compute_stress_tangent(Δϵ, material, state[q_point])
                        # σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, ϕ, η, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d = 
                        state[q_point] = update(material, t, state[q_point], Δϵ, p; kwargs...)
                        # @show state[q_point].ϵ̲̲
                        ϵ = state[q_point].ϵ̲̲
                        σ = state[q_point].σ̲̲
                        # σ_prev, ϵ_prev, ϵᵖ_prev = state[q_point].σ, state[q_point].ϵ, state[q_point].ϵᵖ
                        # ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
                        # σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state[q_point])
                        # Δϵ = ϵ - ϵ_prev
                        dΩ = getdetJdV(cellvalues, q_point)
                        # stateenergy_total[q_point]      = (transpose(ϵ) ⊡ σ) * dΩ
                        # stateenergy_elastic[q_point]    = (transpose(ϵ - state[q_point].ϵᵖ) ⊡ σ) * dΩ
                        # stateenergy_plastic[q_point]    = (transpose(state[q_point].ϵᵖ) ⊡ σ) * dΩ
                        # statetraction[q_point]          = σ # - σ_prev
                        stateenergy_total[q_point]     += (transpose(Δϵ) ⊡ (σ - σ_prev)) * dΩ
                        stateenergy_elastic[q_point]   += (transpose((ϵ - state[q_point].ϵ̲̲⁽ᵖ⁾) - (ϵ_prev - ϵᵖ_prev)) ⊡ (σ - σ_prev)) * dΩ
                        stateenergy_plastic[q_point]   += (transpose(state[q_point].ϵ̲̲⁽ᵖ⁾ - ϵᵖ_prev) ⊡ (σ - σ_prev)) * dΩ
                        statetraction[q_point]         += σ # - σ_prev
                        for i in 1:cell_basefuncs
                            δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
                            # re[i] += (δϵ ⊡ σ) * dΩ # TODO this needs to be uncommented! # add internal force to residual
                            re[i] += ((δϵ - strain[q_point][i]) ⊡ (σ - σ_prev)) * dΩ
                            for j in 1:i # loop only over lower half
                                Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                                # # ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
                                # ke[i, j] += (δϵ - strain[q_point][i]) ⊡ D ⊡ (Δϵ - strain[q_point][j]) * dΩ
                                ke[i, j] += (δϵ - strain[q_point][i]) ⊡ ψ.D⁽ᵉ⁾ ⊡ (Δϵ - strain[q_point][j]) * dΩ
                            end
                            strain[q_point][i] = δϵ
                        end
                    end
                    symmetrize_lower!(ke)
                    assemble!(assembler, eldofs, ke, re)
                end

                # Residual contribution from the Neumann boundary (surface integral)
                # doassemble_neumann!(r, dh, getfacetset(grid, "front_xz"), facetvalues, traction)
                # n_basefuncs = getnbasefunctions(facetvalues)
                rf = zeros(facet_basefuncs)                      # element residual vector
                for fc in FacetIterator(dh, getfacetset(grid, "front_xz"))
                    # Add traction as a negative contribution to the element residual `re`:
                    Ferrite.reinit!(facetvalues, fc)
                    fill!(rf, 0)
                    # rf = @view r[celldofs(fc)]
                    # rf .= 0.0
                    # uf = @view u[celldofs(fc)]
                    stateenergy_total = state_energy_total[:, fc.cc.cellid]
                    stateenergy_total_itr = state_energy_total_itr[:, fc.cc.cellid]
                    stateenergy_total_old = state_energy_total_old[:, fc.cc.cellid]
                    statetraction = state_traction[:, fc.cc.cellid]
                    for q_point in 1:nqp_facet
                        traction_u = (stateenergy_total[q_point] - stateenergy_total_itr[q_point]) * getnormal(facetvalues, q_point)
                        # traction_u = stateenergy_total[q_point] * getnormal(facetvalues, q_point)
                        # traction_u = stateenergy_total[q_point] * Ferrite.Vec{3}([1.0, 0.0, 0.0])
                        # traction_u = (stateenergy_total[q_point] - stateenergy_total_itr[q_point]) * Ferrite.Vec{3}([0.0, 1.0, 0.0])
                        # traction_u = (stateenergy_total[q_point] - stateenergy_total_old[q_point]) * Ferrite.Vec{3}([1.0, 0.0, 0.0])
                        # traction_u = statetraction[q_point] * getnormal(facetvalues, q_point)
                        # traction_u = statetraction[q_point] * Ferrite.Vec{3}([0.0, 1.0, 0.0])
                        # traction_u = Ferrite.Vec{3}(3000e3 .* [1.0, 0.0, 0.0])
                        dΓ = getdetJdV(facetvalues, q_point)
                        # rf[q_point] -= (uf[3(q_point - 1) .+ [1, 2, 3]] ⋅ traction_u) * dΓ
                        for i in 1:facet_basefuncs
                            δu = shape_value(facetvalues, q_point, i)
                            rf[i] -= (δu ⋅ traction_u) * dΓ
                        end
                    end
                    assemble!(r, celldofs(fc), rf)
                end


                state_strain_itr       .= state_strain
                state_energy_total_itr .= state_energy_total
            end

            # break if within tolerance
            norm_r = norm(r[Ferrite.free_dofs(ch)])
            # @printf("\tIteration: %02d \tresidual: %.9f\n", newton_itr, norm_r) # \titx coords: $intersecting_coords") # , $(@sprintf("%.9f", norm_s))")
            if norm_r < NEWTON_TOL
                break
            end


            # increment residual to next NR iteration
            apply_zero!(K, r, ch)
            Δu = Symmetric(K) \ r
            apply_zero!(Δu, ch)
            u -= Δu
            newton_itr += 1
        end

        if newton_itr > NEWTON_M
            error("Reached maximum Newton iterations, aborting")
        end

        # write to results file
            # velocity = zeros(ncells_grid)
            strainᵗ_vonMises    = zeros(gridncells)
            strainᵗ_11          = zeros(gridncells)
            strainᵗ_22          = zeros(gridncells)
            strainᵗ_33          = zeros(gridncells)
            strainᵗ_12          = zeros(gridncells)
            strainᵗ_13          = zeros(gridncells)
            strainᵗ_21          = zeros(gridncells)
            strainᵖ_vonMises    = zeros(gridncells)
            strainᵖ_11          = zeros(gridncells)
            strainᵖ_22          = zeros(gridncells)
            strainᵖ_33          = zeros(gridncells)
            strainᵖ_12          = zeros(gridncells)
            strainᵖ_13          = zeros(gridncells)
            strainᵖ_21          = zeros(gridncells)
            strainᵉ_vonMises    = zeros(gridncells)
            strainᵉ_11          = zeros(gridncells)
            strainᵉ_22          = zeros(gridncells)
            strainᵉ_33          = zeros(gridncells)
            strainᵉ_12          = zeros(gridncells)
            strainᵉ_13          = zeros(gridncells)
            strainᵉ_21          = zeros(gridncells)
            strainrateᵗ_vonMises= zeros(gridncells)
            strainrateᵗ_11      = zeros(gridncells)
            strainrateᵗ_22      = zeros(gridncells)
            strainrateᵗ_33      = zeros(gridncells)
            strainrateᵗ_12      = zeros(gridncells)
            strainrateᵗ_13      = zeros(gridncells)
            strainrateᵗ_21      = zeros(gridncells)
            strainrateᵖ_vonMises= zeros(gridncells)
            strainrateᵖ_11      = zeros(gridncells)
            strainrateᵖ_22      = zeros(gridncells)
            strainrateᵖ_33      = zeros(gridncells)
            strainrateᵖ_12      = zeros(gridncells)
            strainrateᵖ_13      = zeros(gridncells)
            strainrateᵖ_21      = zeros(gridncells)
            strainrateᵉ_vonMises= zeros(gridncells)
            strainrateᵉ_11      = zeros(gridncells)
            strainrateᵉ_22      = zeros(gridncells)
            strainrateᵉ_33      = zeros(gridncells)
            strainrateᵉ_12      = zeros(gridncells)
            strainrateᵉ_13      = zeros(gridncells)
            strainrateᵉ_21      = zeros(gridncells)
            sigma_vonMises      = zeros(gridncells)
            sigma_11            = zeros(gridncells)
            sigma_22            = zeros(gridncells)
            sigma_33            = zeros(gridncells)
            sigma_12            = zeros(gridncells)
            sigma_13            = zeros(gridncells)
            sigma_21            = zeros(gridncells)
            # α_max[i]        = maximum(abs, alpha_vonMises) # hcat(α⃗...)
            # κ_max[i]        = maximum(abs, kappa) # hcat(κ_vec...)
            # ϕ_max[i]        = maximum(abs, phi) # hcat(ϕ⃗...)
            # X_max[i]        = maximum(abs, X) # hcat(X⃗...)
            # d_max[i]        = maximum(abs, d) # hcat(d⃗...)
            alpha_vonMises      = zeros(gridncells)
            alpha_11            = zeros(gridncells)
            alpha_22            = zeros(gridncells)
            alpha_33            = zeros(gridncells)
            alpha_12            = zeros(gridncells)
            alpha_13            = zeros(gridncells)
            alpha_21            = zeros(gridncells)
            kappa               = zeros(gridncells)
            phi                 = zeros(gridncells)
            X                   = zeros(gridncells)
            d                   = zeros(gridncells)
            # κ_values = zeros(gridncells)
            # for (el, state_cells) in enumerate(eachcol(state_material))
            for ((el, state_cells), state_old_cells) in zip(enumerate(eachcol(state_material)), eachcol(state_material_old))
                for (state, state_old) in zip(state_cells, state_old_cells)
                    # velocity[el] += first(rand(n_timesteps))
                    strainᵗ_vonMises[el]    += vonMises(state.ϵ̲̲)
                    strainᵗ_11[el]          += state.ϵ̲̲[1, 1]
                    strainᵗ_22[el]          += state.ϵ̲̲[2, 2]
                    strainᵗ_33[el]          += state.ϵ̲̲[3, 3]
                    strainᵗ_12[el]          += state.ϵ̲̲[1, 2]
                    strainᵗ_13[el]          += state.ϵ̲̲[1, 3]
                    strainᵗ_21[el]          += state.ϵ̲̲[2, 1]
                    strainᵖ_vonMises[el]    += vonMises(state.ϵ̲̲⁽ᵖ⁾)
                    strainᵖ_11[el]          += state.ϵ̲̲⁽ᵖ⁾[1, 1]
                    strainᵖ_22[el]          += state.ϵ̲̲⁽ᵖ⁾[2, 2]
                    strainᵖ_33[el]          += state.ϵ̲̲⁽ᵖ⁾[3, 3]
                    strainᵖ_12[el]          += state.ϵ̲̲⁽ᵖ⁾[1, 2]
                    strainᵖ_13[el]          += state.ϵ̲̲⁽ᵖ⁾[1, 3]
                    strainᵖ_21[el]          += state.ϵ̲̲⁽ᵖ⁾[2, 1]
                    strainᵉ_vonMises[el]    += (vonMises(state.ϵ̲̲) - vonMises(state.ϵ̲̲⁽ᵖ⁾))
                    strainᵉ_11[el]          += (state.ϵ̲̲[1, 1] - state.ϵ̲̲⁽ᵖ⁾[1, 1])
                    strainᵉ_22[el]          += (state.ϵ̲̲[2, 2] - state.ϵ̲̲⁽ᵖ⁾[2, 2])
                    strainᵉ_33[el]          += (state.ϵ̲̲[3, 3] - state.ϵ̲̲⁽ᵖ⁾[3, 3])
                    strainᵉ_12[el]          += (state.ϵ̲̲[1, 2] - state.ϵ̲̲⁽ᵖ⁾[1, 2])
                    strainᵉ_13[el]          += (state.ϵ̲̲[1, 3] - state.ϵ̲̲⁽ᵖ⁾[1, 3])
                    strainᵉ_21[el]          += (state.ϵ̲̲[2, 1] - state.ϵ̲̲⁽ᵖ⁾[2, 1])
                    strainrateᵗ_vonMises[el]+= (vonMises(state.ϵ̲̲) - vonMises(state_old.ϵ̲̲)) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵗ_11[el]      += (state.ϵ̲̲[1, 1] - state_old.ϵ̲̲[1, 1]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵗ_22[el]      += (state.ϵ̲̲[2, 2] - state_old.ϵ̲̲[2, 2]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵗ_33[el]      += (state.ϵ̲̲[3, 3] - state_old.ϵ̲̲[3, 3]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵗ_12[el]      += (state.ϵ̲̲[1, 2] - state_old.ϵ̲̲[1, 2]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵗ_13[el]      += (state.ϵ̲̲[1, 3] - state_old.ϵ̲̲[1, 3]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵗ_21[el]      += (state.ϵ̲̲[2, 1] - state_old.ϵ̲̲[2, 1]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵖ_vonMises[el]+= (vonMises(state.ϵ̲̲⁽ᵖ⁾) - vonMises(state_old.ϵ̲̲⁽ᵖ⁾)) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵖ_11[el]      += (state.ϵ̲̲⁽ᵖ⁾[1, 1] - state_old.ϵ̲̲[1, 1]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵖ_22[el]      += (state.ϵ̲̲⁽ᵖ⁾[2, 2] - state_old.ϵ̲̲[2, 2]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵖ_33[el]      += (state.ϵ̲̲⁽ᵖ⁾[3, 3] - state_old.ϵ̲̲[3, 3]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵖ_12[el]      += (state.ϵ̲̲⁽ᵖ⁾[1, 2] - state_old.ϵ̲̲[1, 2]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵖ_13[el]      += (state.ϵ̲̲⁽ᵖ⁾[1, 3] - state_old.ϵ̲̲[1, 3]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵖ_21[el]      += (state.ϵ̲̲⁽ᵖ⁾[2, 1] - state_old.ϵ̲̲[2, 1]) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵉ_vonMises[el]+= ((vonMises(state.ϵ̲̲) - vonMises(state.ϵ̲̲⁽ᵖ⁾)) - (vonMises(state_old.ϵ̲̲) - vonMises(state_old.ϵ̲̲⁽ᵖ⁾))) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵉ_11[el]      += ((state.ϵ̲̲[1, 1] - state.ϵ̲̲⁽ᵖ⁾[1, 1]) - (state_old.ϵ̲̲[1, 1] - state_old.ϵ̲̲⁽ᵖ⁾[1, 1])) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵉ_22[el]      += ((state.ϵ̲̲[2, 2] - state.ϵ̲̲⁽ᵖ⁾[2, 2]) - (state_old.ϵ̲̲[2, 2] - state_old.ϵ̲̲⁽ᵖ⁾[2, 2])) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵉ_33[el]      += ((state.ϵ̲̲[3, 3] - state.ϵ̲̲⁽ᵖ⁾[3, 3]) - (state_old.ϵ̲̲[3, 3] - state_old.ϵ̲̲⁽ᵖ⁾[3, 3])) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵉ_12[el]      += ((state.ϵ̲̲[1, 2] - state.ϵ̲̲⁽ᵖ⁾[1, 2]) - (state_old.ϵ̲̲[1, 2] - state_old.ϵ̲̲⁽ᵖ⁾[1, 2])) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵉ_13[el]      += ((state.ϵ̲̲[1, 3] - state.ϵ̲̲⁽ᵖ⁾[1, 3]) - (state_old.ϵ̲̲[1, 3] - state_old.ϵ̲̲⁽ᵖ⁾[1, 3])) / ψ.Δt # (dN * ψ.Δt)
                    strainrateᵉ_21[el]      += ((state.ϵ̲̲[2, 1] - state.ϵ̲̲⁽ᵖ⁾[2, 1]) - (state_old.ϵ̲̲[2, 1] - state_old.ϵ̲̲⁽ᵖ⁾[2, 1])) / ψ.Δt # (dN * ψ.Δt)
                    # @show state.σ̲̲
                    sigma_vonMises[el]      += vonMises(state.σ̲̲)
                    sigma_11[el]            += state.σ̲̲[1, 1]
                    sigma_22[el]            += state.σ̲̲[2, 2]
                    sigma_33[el]            += state.σ̲̲[3, 3]
                    sigma_12[el]            += state.σ̲̲[1, 2]
                    sigma_13[el]            += state.σ̲̲[1, 3]
                    sigma_21[el]            += state.σ̲̲[2, 1]
                    alpha_vonMises[el]      += vonMises(state.α̲̲)
                    alpha_11[el]            += state.α̲̲[1, 1]
                    alpha_22[el]            += state.α̲̲[2, 2]
                    alpha_33[el]            += state.α̲̲[3, 3]
                    alpha_12[el]            += state.α̲̲[1, 2]
                    alpha_13[el]            += state.α̲̲[1, 3]
                    alpha_21[el]            += state.α̲̲[2, 1]
                    kappa[el]               += state.κ
                    phi[el]                 += state.ϕ
                    X[el]                   += state.X
                    d[el]                   += state.d
                    # κ_values[el] += state.k*material_endplate.H
                end
                # velocity[el] /= length(cell_states)
                strainᵗ_vonMises[el]    /= length(state_cells)
                strainᵗ_11[el]          /= length(state_cells)
                strainᵗ_22[el]          /= length(state_cells)
                strainᵗ_33[el]          /= length(state_cells)
                strainᵗ_12[el]          /= length(state_cells)
                strainᵗ_13[el]          /= length(state_cells)
                strainᵗ_21[el]          /= length(state_cells)
                strainᵖ_vonMises[el]    /= length(state_cells)
                strainᵖ_11[el]          /= length(state_cells)
                strainᵖ_22[el]          /= length(state_cells)
                strainᵖ_33[el]          /= length(state_cells)
                strainᵖ_12[el]          /= length(state_cells)
                strainᵖ_13[el]          /= length(state_cells)
                strainᵖ_21[el]          /= length(state_cells)
                strainᵉ_vonMises[el]    /= length(state_cells)
                strainᵉ_11[el]          /= length(state_cells)
                strainᵉ_22[el]          /= length(state_cells)
                strainᵉ_33[el]          /= length(state_cells)
                strainᵉ_12[el]          /= length(state_cells)
                strainᵉ_13[el]          /= length(state_cells)
                strainᵉ_21[el]          /= length(state_cells)
                strainrateᵗ_vonMises[el]/= length(state_cells)
                strainrateᵗ_11[el]      /= length(state_cells)
                strainrateᵗ_22[el]      /= length(state_cells)
                strainrateᵗ_33[el]      /= length(state_cells)
                strainrateᵗ_12[el]      /= length(state_cells)
                strainrateᵗ_13[el]      /= length(state_cells)
                strainrateᵗ_21[el]      /= length(state_cells)
                strainrateᵖ_vonMises[el]/= length(state_cells)
                strainrateᵖ_11[el]      /= length(state_cells)
                strainrateᵖ_22[el]      /= length(state_cells)
                strainrateᵖ_33[el]      /= length(state_cells)
                strainrateᵖ_12[el]      /= length(state_cells)
                strainrateᵖ_13[el]      /= length(state_cells)
                strainrateᵖ_21[el]      /= length(state_cells)
                strainrateᵉ_vonMises[el]/= length(state_cells)
                strainrateᵉ_11[el]      /= length(state_cells)
                strainrateᵉ_22[el]      /= length(state_cells)
                strainrateᵉ_33[el]      /= length(state_cells)
                strainrateᵉ_12[el]      /= length(state_cells)
                strainrateᵉ_13[el]      /= length(state_cells)
                strainrateᵉ_21[el]      /= length(state_cells)
                sigma_vonMises[el]      /= length(state_cells)
                sigma_11[el]            /= length(state_cells)
                sigma_22[el]            /= length(state_cells)
                sigma_33[el]            /= length(state_cells)
                sigma_12[el]            /= length(state_cells)
                sigma_13[el]            /= length(state_cells)
                sigma_21[el]            /= length(state_cells)
                alpha_vonMises[el]      /= length(state_cells)
                alpha_11[el]            /= length(state_cells)
                alpha_22[el]            /= length(state_cells)
                alpha_33[el]            /= length(state_cells)
                alpha_12[el]            /= length(state_cells)
                alpha_13[el]            /= length(state_cells)
                alpha_21[el]            /= length(state_cells)
                kappa[el]               /= length(state_cells)
                phi[el]                 /= length(state_cells)
                X[el]                   /= length(state_cells)
                d[el]                   /= length(state_cells)
                # κ_values[el] /= length(state_cells)
            end
            strain_max[i] = maximum(abs, strainᵗ_vonMises) # maximum displacement in current timestep
            strainrate_max[i] = maximum(abs, strainrateᵗ_vonMises) # maximum displacement in current timestep
            stress_max[i] = maximum(abs, sigma_vonMises) # maximum displacement in current timestep
            α_max[i]        = maximum(abs, alpha_vonMises) # hcat(α⃗...)
            κ_max[i]        = maximum(abs, kappa) # hcat(κ_vec...)
            ϕ_max[i]        = maximum(abs, phi) # hcat(ϕ⃗...)
            X_max[i]        = maximum(abs, X) # hcat(X⃗...)
            d_max[i]        = maximum(abs, d) # hcat(d⃗...)
            # strain_max[timestep] = maximum(abs, strainᵗ_22) # maximum displacement in current timestep
            # stress_max[timestep] = maximum(abs, sigma_22) # maximum displacement in current timestep
            Vₐ = try
                derive_serial(collect(range(0, 0.1, M))[begin:Int64((0.02 / (0.1 / (M - 1))) ÷ 1 + 1)][begin:i],
                    ΔX, 1:i, i, :five, 1, ψ.Δt)
            catch
                try
                    (ΔX[i] - ΔX[i - 1]) / ψ.Δt
                catch
                    V
                end
            end
            velocity_applied = zeros((3, gridnnodes))
            velocity_actual = zeros((3, gridnnodes))
            for i in getnodeset(grid, "front_xz")
                velocity_applied[:, i] = V .* [1.0, 1.0, 0.0]
                velocity_actual[:, i] = Vₐ .* [1.0, 1.0, 0.0]
            end
            # velocity = zeros(n_dofs)
            # for i in getnodeset(grid, "indenter_tip")
            #     velocity[vectorstrafe_i(i)] .= V .* [0., 0., 1.]
            # end
            VTKGridFile("rve" * "-t$i", dh) do vtk
                # vtk["TimeValue"] = float(timestep)
                write_solution(vtk, dh, u) # displacement field
                # // DONE [20250110] (JMA3): figure out how to edit output vectors for visualization
                # [20250206T0911] (JMA3): `write_node_data` ≡ `write_solution` for display vectors in ParaView
                # write_solution(vtk, dh, fill(t, n_dofs), "v") # velocity field
                # write_cell_data(vtk, velocity, "v") # velocity field
                write_node_data(vtk, velocity_applied, "Velocity (Applied)")
                write_node_data(vtk, velocity_actual, "Velocity (Actual)")
                # write_solution(vtk, dh, u) # displacement field

                write_cell_data(vtk, strainᵗ_vonMises,      "strain_t")
                write_cell_data(vtk, strainᵗ_11,            "strain_t_11")
                write_cell_data(vtk, strainᵗ_22,            "strain_t_22")
                write_cell_data(vtk, strainᵗ_33,            "strain_t_33")
                write_cell_data(vtk, strainᵗ_12,            "strain_t_12")
                write_cell_data(vtk, strainᵗ_13,            "strain_t_13")
                write_cell_data(vtk, strainᵗ_21,            "strain_t_21")
                write_cell_data(vtk, strainᵖ_vonMises,      "strain_p")
                write_cell_data(vtk, strainᵖ_11,            "strain_p_11")
                write_cell_data(vtk, strainᵖ_22,            "strain_p_22")
                write_cell_data(vtk, strainᵖ_33,            "strain_p_33")
                write_cell_data(vtk, strainᵖ_12,            "strain_p_12")
                write_cell_data(vtk, strainᵖ_13,            "strain_p_13")
                write_cell_data(vtk, strainᵖ_21,            "strain_p_21")
                write_cell_data(vtk, strainᵉ_vonMises,      "strain_e")
                write_cell_data(vtk, strainᵉ_11,            "strain_e_11")
                write_cell_data(vtk, strainᵉ_22,            "strain_e_22")
                write_cell_data(vtk, strainᵉ_33,            "strain_e_33")
                write_cell_data(vtk, strainᵉ_12,            "strain_e_12")
                write_cell_data(vtk, strainᵉ_13,            "strain_e_13")
                write_cell_data(vtk, strainᵉ_21,            "strain_e_21")
                write_cell_data(vtk, strainrateᵗ_vonMises,  "strainrate_t")
                write_cell_data(vtk, strainrateᵗ_11,        "strainrate_t_11")
                write_cell_data(vtk, strainrateᵗ_22,        "strainrate_t_22")
                write_cell_data(vtk, strainrateᵗ_33,        "strainrate_t_33")
                write_cell_data(vtk, strainrateᵗ_12,        "strainrate_t_12")
                write_cell_data(vtk, strainrateᵗ_13,        "strainrate_t_13")
                write_cell_data(vtk, strainrateᵗ_21,        "strainrate_t_21")
                write_cell_data(vtk, strainrateᵖ_vonMises,  "strainrate_p")
                write_cell_data(vtk, strainrateᵖ_11,        "strainrate_p_11")
                write_cell_data(vtk, strainrateᵖ_22,        "strainrate_p_22")
                write_cell_data(vtk, strainrateᵖ_33,        "strainrate_p_33")
                write_cell_data(vtk, strainrateᵖ_12,        "strainrate_p_12")
                write_cell_data(vtk, strainrateᵖ_13,        "strainrate_p_13")
                write_cell_data(vtk, strainrateᵖ_21,        "strainrate_p_21")
                write_cell_data(vtk, strainrateᵉ_vonMises,  "strainrate_e")
                write_cell_data(vtk, strainrateᵉ_11,        "strainrate_e_11")
                write_cell_data(vtk, strainrateᵉ_22,        "strainrate_e_22")
                write_cell_data(vtk, strainrateᵉ_33,        "strainrate_e_33")
                write_cell_data(vtk, strainrateᵉ_12,        "strainrate_e_12")
                write_cell_data(vtk, strainrateᵉ_13,        "strainrate_e_13")
                write_cell_data(vtk, strainrateᵉ_21,        "strainrate_e_21")
                write_cell_data(vtk, sigma_vonMises,        "sigma_t")
                write_cell_data(vtk, sigma_11,              "sigma_11")
                write_cell_data(vtk, sigma_22,              "sigma_22")
                write_cell_data(vtk, sigma_33,              "sigma_33")
                write_cell_data(vtk, sigma_12,              "sigma_12")
                write_cell_data(vtk, sigma_13,              "sigma_13")
                write_cell_data(vtk, sigma_21,              "sigma_21")
                # for i in 1:3, j in 1:3
                #     σij = [x[i, j] for x in σ]
                #     write_cell_data(vtk, σij, "sigma_$(i)$(j)")
                # end
                # write_cell_data(vtk, κ_values, "Drag stress [Pa]")
                # Ferrite.write_constraints(vtk, ch)
                # pvd[timestep] = vtk
                collection_add_timestep(pvd, vtk, t)
            end # closes vtk
        #

        # update the old states with the converged values for next timestep
        u₀ = u
        V += ΔV
        # states_material_old        .= state_material
        # states_strain_rate         .= (state_strain - state_strain_old) / ψ.Δt # (dN * ψ.Δt)
        # states_strain_old          .= state_strain
        # # states_energy_total_old    .= state_energy_total
        # # states_energy_elastic_old  .= state_energy_elastic
        # # states_energy_plastic_old  .= state_energy_plastic
        # states_energy_total[timestep]   = state_energy_total
        # states_energy_elastic[timestep] = state_energy_elastic
        # states_energy_plastic[timestep] = state_energy_plastic
        states_material[i]       = state_material
        states_strain_rate[i]    = (state_strain - state_strain_old) / ψ.Δt # (dN * ψ.Δt)
        states_energy_total[i]   = state_energy_total
        states_energy_elastic[i] = state_energy_elastic
        states_energy_plastic[i] = state_energy_plastic
        u_max[i] = maximum(abs, u) # maximum displacement in current timestep
        ###########################################################
        # # println("")
        # # # @show i, t, ψ.θ, ψ.Δt, d
        # σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, κₛ, ϕ, η, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d = update(ψ, t, σ̲̲, ϵ̲̲, ϵ̲̲⁽ᵖ⁾, α̲̲, κ, κₛ, Si, ϕ, damirr, η, νᵥ, ϕ̇, X, XR, XH, Xd, Xs, d, p; kwargs...)
        # push!(ϵ⃗, ϵ̲̲)
        # push!(σ⃗, σ̲̲)
        # push!(α⃗, sqrt_threehalves * norm_symvec(α̲̲))
        # push!(κ_vec, κ)
        # push!(ϕ⃗, ϕ)
        # push!(X⃗, X)
        # push!(d⃗, d)
    end
    # return (data=(
    #     ϵ=hcat(ϵ⃗...), σ=hcat(σ⃗...),
    #     α=hcat(α⃗...), κ=hcat(κ_vec...),
    #     ϕ=hcat(ϕ⃗...), X=hcat(X⃗...),
    #     d=hcat(d⃗...),
    # ),)
    return (data=(
        ϵ=hcat(strain_max...),
        σ=hcat(stress_max...),
        α=hcat(α_max...), # hcat(α⃗...)
        κ=hcat(κ_max...), # hcat(κ_vec...)
        ϕ=hcat(ϕ_max...), # hcat(ϕ⃗...)
        X=hcat(X_max...), # hcat(X⃗...)
        d=hcat(d_max...), # hcat(d⃗...)
        # α=hcat(α⃗...), κ=hcat(κ_vec...),
        # ϕ=hcat(ϕ⃗...), X=hcat(X⃗...),
        # d=hcat(d⃗...),
    ),)
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
ContinuumMechanicsBase.parameters(::Cho2019UnifiedStaticDynamicTensor) = (
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
                        :NK,    # * [20250402T1521] (JMA3): I think this is the modifier for finding the k-root
                                # *                         (see Eq. 4.22 in HEC dissertation)
                                # *                         (c. f. `optimize.py` that NK=2.0 by default)
                                # ! [20250422T1121] (JMA3): This is exponent on κ in rate equation.
                                # !                         Bammann assumed 2 for metals for dislocation creep in Power Law
    ## torsion, tension/compression
    :C₁₉, :C₂₀,
    ## dynamic recrystallization
    :Cx1, :Cx2, :Cdp,
    :Cx3, :Cx4, :Csp,
    :Cx5, :Cxa, :Cxb, :Cxc,
    ## static RX (grain growth)
    # :n, :ω₀, # E⁺, V⁺, R,
    ## grain size
    # d₀, Cg1, Cg2, Cg3, z,
    :Cg1, :Cg2, :Cg3, # :z,
    # :Cg1, :Cg2, :Cg3, :Cg4, # :z,
    ## damage
    ### nucleation
    # 𝒹, 𝒻, Kic, a, b, c,
    :Dc, :a, :b, :c,
    # Cnuc, Tnuc, R₀, nn, Tgrw,
    :Cnuc, :Tnuc, :nn, :Tgrw,
    ## irradiation hardening
    :kr1, :krt, :kr2, :kr3, :kp1, :kpt, :kp2
)

nothing
# end # end of module
