abstract type BCJ end
abstract type BCJMetal              <: BCJ end
abstract type Bammann1990Modeling   <: BCJMetal end
abstract type DK                    <: BCJMetal end
abstract type ISVMetal{T<:BCJMetal} end
abstract type KinematicHardening{T} <: ISVMetal{T} end # α__
abstract type IsotropicHardening{T} <: ISVMetal{T} end # κ
abstract type Damage{T}             <: ISVMetal{T} end # ϕ

symmetricmagnitude(tensor::Vector{<:Real}) = √( sum(tensor[1:3] .^ 2.) + 2sum(tensor[4:6] .^ 2.) )
function symmetricvonMises(tensor::Union{Vector{<:Real}, SubArray{<:Real}})::AbstractFloat
    σvM = sum(map(x->x^2., [tensor[1] - tensor[2], tensor[2] - tensor[3], tensor[3] - tensor[1]])) + (
        6sum(map(x->x^2., [tensor[4], tensor[5], tensor[6]])))
    return √(σvM / 2.)
end
function symmetricvonMises(tensor::Matrix{<:Real})::Vector{AbstractFloat}
    return map(symmetricvonMises, eachcol(tensor))
end

struct BCJMetalStrainControl{T1<:Integer, T2<:AbstractFloat} <: BCJMetal
    θ       ::T2                # applied temperature
    ϵ_dot   ::T2                # applied strain rate
    ϵₙ      ::T2                # final strain
    N       ::T1                # number of strain increments
    istate  ::T1                # load type (1: uniaxial tension; 2: torsion)
    params  ::Dict{String, T2}  # material constants
end

mutable struct BCJMetalCurrentConfiguration{Version<:BCJMetal, T<:AbstractFloat} <: BCJMetal
    N               ::Integer   # number of strain increments
    θ               ::T         # applied temperature
    μ               ::T         # shear modulus at temperature, θ
    σ__             ::Vector{T} # deviatoric stress tensor
    ϵₚ__            ::Vector{T} # plastic strain tensor
    ϵ_dot_plastic__ ::Vector{T} # plastic strain rate
    ϵ__             ::Vector{T} # total strain tensor
    ϵ_dot_effective ::T         # strain rate (effective)
    Δϵ              ::Vector{T} # total strain tensor step
    Δt              ::T         # time step
    V               ::T         # strain rate sensitivity of yield stress at temperature, θ
    Y               ::T         # rate independent yield stress at temperature, θ
    f               ::T         # strain rate at which yield becomes strain rate dependent at temperature, θ
    h               ::T         # kinematic hardening modulus at temperature, θ
    r_d             ::T         # dynamic recovery of kinematic hardening at temperature, θ
    r_s             ::T         # diffusion controlled static/thermal recovery of kinematic hardening at temperature, θ
    H               ::T         # isotropic hardening modulus at temperature, θ
    R_d             ::T         # dynamic recovery of isotropic hardening at temperature, θ
    R_s             ::T         # diffusion controlled static/thermal recovery of isotropic hardening at temperature, θ
    α__             ::Vector{T} # kinematic hardening tensor
    κ               ::T         # isotropic hardening scalar
    β               ::T         # yield function
    ξ__             ::Vector{T} # overstress tensor (S - 2/3*alpha)
    σₜᵣ__           ::Vector{T} # deviatoric stress tensor (trial)
    αₜᵣ__           ::Vector{T} # kinematic hardening tensor (trial)
    κₜᵣ             ::T         # isotropic hardening (trial)
end

mutable struct BCJMetalConfigurationHistory{T<:AbstractFloat} <: BCJMetal
    σ__             ::Matrix{T} # deviatoric stress tensor
    ϵₚ__            ::Matrix{T} # plastic strain tensor
    ϵ_dot_plastic__ ::Matrix{T} # plastic strain rate
    ϵ__             ::Matrix{T} # total strain tensor
    α__             ::Matrix{T} # kinematic hardening tensor
    κ               ::Vector{T} # isotropic hardening scalar
    ξ__             ::Matrix{T} # overstress tensor (S - 2/3*alpha)
end

function Base.:+(x::T, y::T) where {T<:BCJMetalConfigurationHistory}
    return BCJMetalConfigurationHistory{eltype(x.σ__)}(
        hcat(x.σ__,                y.σ__),
        hcat(x.ϵₚ__,               y.ϵₚ__),
        hcat(x.ϵ_dot_plastic__,    y.ϵ_dot_plastic__),
        hcat(x.ϵ__,                y.ϵ__ .+ x.ϵ__[:, end]),
        hcat(x.α__,                y.α__),
        vcat(x.κ,                  y.κ),
        hcat(x.ξ__,                y.ξ__)
    )
end

function Base.copyto!(reference::BCJMetalCurrentConfiguration, history::BCJMetalConfigurationHistory)
    # for attr ∈ (:θ, :V, :Y, :f, :h, :r_d, :r_s, :H, :R_d, :R_s, :α__, :κ, :β, :ξ__)
    #     setfield!(reference, attr, getfield(current, attr))
    # end
    reference.σ__               = history.σ__[:, end]
    # reference.ϵₚ__              = history.ϵₚ__[:, end]
    # reference.ϵ_dot_plastic__   = history.ϵ_dot_plastic__[:, end]
    # reference.ϵ__               = history.ϵ__[:, end]
    reference.α__               = history.α__[:, end]
    reference.κ                 = history.κ[end]
    reference.ξ__               = history.ξ__[:, end]
    return nothing
end

function record!(history::BCJMetalConfigurationHistory, i::Integer, current::BCJMetalCurrentConfiguration)
    history.σ__[:, i]              .= current.σ__
    history.ϵₚ__[:, i]             .= current.ϵₚ__
    history.ϵ_dot_plastic__[:, i]  .= current.ϵ_dot_plastic__
    history.ϵ__[:, i]              .= current.ϵ__
    history.α__[:, i]              .= current.α__
    history.κ[i]                    = current.κ
    history.ξ__[:, i]              .= current.ξ__
    return nothing
end

function bcjmetalreferenceconfiguration(::Type{Bammann1990Modeling}, bcj::BCJMetalStrainControl)::Tuple{BCJMetalCurrentConfiguration, BCJMetalCurrentConfiguration, BCJMetalConfigurationHistory}
    θ       = bcj.θ
    ϵ_dot   = bcj.ϵ_dot
    ϵₙ      = bcj.ϵₙ
    N       = bcj.N
    istate  = bcj.istate
    params  = bcj.params
    M       = N + 1
    T       = typeof(float(θ))
    # breakout params into easy variables
    G       = params["bulk_mod"]
    μ       = params["shear_mod"]
    # params_keys = keys(params)
    # C = params[params_keys[findall(r"C\d+", params_keys)]]
    C1      = params["C01"]
    C2      = params["C02"]
    C3      = params["C03"]
    C4      = params["C04"]
    C5      = params["C05"]
    C6      = params["C06"]
    C7      = params["C07"]
    C8      = params["C08"]
    C9      = params["C09"]
    C10     = params["C10"]
    C11     = params["C11"]
    C12     = params["C12"]
    C13     = params["C13"]
    C14     = params["C14"]
    C15     = params["C15"]
    C16     = params["C16"]
    C17     = params["C17"]
    C18     = params["C18"]
    C19     = params["C19"]
    C20     = params["C20"]


    # array declarations
    # * tenXirs: # = [#_11, #_22, #_33, #_12, #_23, #_13]
    ## OSVs
    σ__             = zeros(T, 6)   # deviatoric stress
    ϵₚ__            = zeros(T, 6)   # plastic strain
    ϵ_dot_plastic__ = zeros(T, 6)   # plastic strain rate
    ϵ__             = zeros(T, 6)   # total strain
    ## ISVs
    α__             = fill(1e-7, 6) # alpha: kinematic hardening
    κ               = 0.            # kappa: isotropic hardening
    ## holding values
    Δϵ              = zeros(T, 6)   # strain increment
    σₜᵣ__           = zeros(T, 6)   # trial stress  (deviatoric)
    αₜᵣ__           = zeros(T, 6)   # trial kinematic
    κₜᵣ             = 0.
    ξ__             = zeros(T, 6)   # overstress (S - 2/3*alpha)

    # state evaluation - loading type
    ϵ_dot_effective = if istate == 1    # uniaxial tension
        δϵ  = ϵₙ / N
        Δϵ .= [δϵ, -0.499δϵ, -0.499δϵ, 0., 0., 0.]
        Δt  = δϵ / ϵ_dot # timestep
        ϵ_dot
    elseif istate == 2                  # torsion
        # convert equivalent strain to true shear strain
        ϵₙ *= 0.5 * √(3.)
        Δϵ .= [0., 0., 0., ϵₙ / N, 0., 0.]
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[3] / ϵ_dot            # timestep
        2ϵ_dot / √3.
    end


    # temperature dependent constants
    V   = C1    * exp( -C2 / θ )
    Y   = C3    * exp(  C4 / θ )
    f   = C5    * exp( -C6 / θ )

    β   = Y + (V * asinh( ϵ_dot_effective / f ))

    r_d = C7    * exp( -C8  / θ )
    h   = C9    * exp(  C10 * θ )
    r_s = C11   * exp( -C12 / θ )

    R_d = C13   * exp( -C14 / θ )
    H   = C15   * exp(  C16 * θ )
    R_s = C17   * exp( -C18 / θ )
    current = BCJMetalCurrentConfiguration{Bammann1990Modeling, T}(N, θ, μ,
        σ__, ϵₚ__, ϵ_dot_plastic__, ϵ__, ϵ_dot_effective, Δϵ, Δt,
        V, Y, f, h, r_d, r_s, H, R_d, R_s, α__, κ, β, ξ__, σₜᵣ__, αₜᵣ__, κₜᵣ)
    history = BCJMetalConfigurationHistory{T}(
        ## OSVs
        Matrix{T}(undef, (6, M)),   # deviatoric stress, σ__
        Matrix{T}(undef, (6, M)),   # plastic strain, ϵₚ__
        Matrix{T}(undef, (6, M)),   # plastic strain rate, ϵ_dot_plastic__
        Matrix{T}(undef, (6, M)),   # total strain, ϵ__
        ## ISVs
        Matrix{T}(undef, (6, M)),   # kinematic hardening, α__
        Vector{T}(undef,     M ),   # isotropic hardening, κ
        ## holding values
        Matrix{T}(undef, (6, M))    # overstress (S - 2/3*alpha), ξ__
    )
    record!(history, 1, current)
    return (current, current, history)
end

function bcjmetalreferenceconfiguration(::Type{DK}, bcj::BCJMetalStrainControl)::Tuple{BCJMetalCurrentConfiguration, BCJMetalCurrentConfiguration, BCJMetalConfigurationHistory}
    θ       = bcj.θ
    ϵ_dot   = bcj.ϵ_dot
    ϵₙ      = bcj.ϵₙ
    N       = bcj.N
    istate  = bcj.istate
    params  = bcj.params
    M       = N + 1
    T       = typeof(float(θ))
    # breakout params into easy variables
    G       = params["bulk_mod"]
    μ       = params["shear_mod"]
    # params_keys = keys(params)
    # C = params[params_keys[findall(r"C\d+", params_keys)]]
    C1      = params["C01"]
    C2      = params["C02"]
    C3      = params["C03"]
    C4      = params["C04"]
    C5      = params["C05"]
    C6      = params["C06"]
    C7      = params["C07"]
    C8      = params["C08"]
    C9      = params["C09"]
    C10     = params["C10"]
    C11     = params["C11"]
    C12     = params["C12"]
    C13     = params["C13"]
    C14     = params["C14"]
    C15     = params["C15"]
    C16     = params["C16"]
    C17     = params["C17"]
    C18     = params["C18"]
    C19     = params["C19"]
    C20     = params["C20"]


    # array declarations
    # * tenXirs: # = [#_11, #_22, #_33, #_12, #_23, #_13]
    ## OSVs
    σ__             = zeros(T, 6)   # deviatoric stress
    ϵₚ__            = zeros(T, 6)   # plastic strain
    ϵ_dot_plastic__ = zeros(T, 6)   # plastic strain rate
    ϵ__             = zeros(T, 6)   # total strain
    ## ISVs
    α__             = fill(1e-7, 6) # alpha: kinematic hardening
    κ               = 0.            # kappa: isotropic hardening
    ## holding values
    Δϵ              = zeros(T, 6)   # strain increment
    σₜᵣ__           = zeros(T, 6)   # trial stress  (deviatoric)
    αₜᵣ__           = zeros(T, 6)   # trial kinematic
    κₜᵣ             = 0.
    ξ__             = zeros(T, 6)   # overstress (S - 2/3*alpha)


    # state evaluation - loading type
    ϵ_dot_effective = if istate == 1    # uniaxial tension
        δϵ  = ϵₙ / N
        Δϵ .= [δϵ, -0.499δϵ, -0.499δϵ, 0., 0., 0.]
        Δt  = δϵ / ϵ_dot # timestep
        ϵ_dot
    elseif istate == 2                  # torsion
        # convert equivalent strain to true shear strain
        ϵₙ *= 0.5 * √(3.)
        Δϵ .= [0., 0., 0., ϵₙ / N, 0., 0.]
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[3] / ϵ_dot            # timestep
        2ϵ_dot / √3.
    end


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
    current = BCJMetalCurrentConfiguration{DK, T}(N, θ, μ,
        σ__, ϵₚ__, ϵ_dot_plastic__, ϵ__, ϵ_dot_effective, Δϵ, Δt,
        V, Y, f, h, r_d, r_s, H, R_d, R_s, α__, κ, β, ξ__, σₜᵣ__, αₜᵣ__, κₜᵣ)
    history = BCJMetalConfigurationHistory{T}(
        ## OSVs
        Matrix{T}(undef, (6, M)),   # deviatoric stress, σ__
        Matrix{T}(undef, (6, M)),   # plastic strain, ϵₚ__
        Matrix{T}(undef, (6, M)),   # plastic strain rate, ϵ_dot_plastic__
        Matrix{T}(undef, (6, M)),   # total strain, ϵ__
        ## ISVs
        Matrix{T}(undef, (6, M)),   # kinematic hardening, α__
        Vector{T}(undef,     M ),   # isotropic hardening, κ
        ## holding values
        Matrix{T}(undef, (6, M))    # overstress (S - 2/3*alpha), ξ__
    )
    record!(history, 1, current)
    return (current, current, history)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function solve!(bcj::BCJMetalCurrentConfiguration{Bammann1990Modeling, <:AbstractFloat},
        history::BCJMetalConfigurationHistory)
    μ, Δϵ, Δt       = bcj.μ, bcj.Δϵ, bcj.Δt
    ϵ_dot_effective = bcj.ϵ_dot_effective
    V, Y, f, β      = bcj.V, bcj.Y, bcj.f, bcj.β
    h, r_d, r_s     = bcj.h, bcj.r_d, bcj.r_s # alpha
    H, R_d, R_s     = bcj.H, bcj.R_d, bcj.R_s # kappa
    # timestep calculations
    for i ∈ range(2, bcj.N + 1)
        α_mag       = symmetricmagnitude(bcj.α__)
        # trial guesses: ISVs (from recovery) and stress
        recovery    = Δt * (r_d * ϵ_dot_effective + r_s) * α_mag    # recovery for alpha (kinematic hardening)
        Recovery    = Δt * (R_d * ϵ_dot_effective + R_s) * bcj.κ    # recovery for kappa (isotropic hardening)
        αₜᵣ__       = bcj.α__  .* (1 - recovery)
        κₜᵣ         = bcj.κ     * (1 - Recovery)

        ## trial stress guess
        σₜᵣ__       = bcj.σ__ + (2μ .* Δϵ)           # trial stress
        bcj.ξ__    .= σₜᵣ__ - αₜᵣ__                 # trial overstress original
        # ξ__          .= σₜᵣ__ - sqrt23 .* αₜᵣ__   # trial overstress FIT
        ξ_mag       = symmetricmagnitude(bcj.ξ__)



        # ----------------------------------- #
        ###   ---   YIELD CRITERION   ---   ###
        # ----------------------------------- #
        flow_rule = ξ_mag - κₜᵣ - β         # same as vumat20
        # Crit = Xi_mag - (Katr + β) #changed to FIT
        if flow_rule <= 0.      # elastic
            # trial guesses are correct
            bcj.σ__    .= σₜᵣ__
            bcj.α__    .= αₜᵣ__
            bcj.κ       = κₜᵣ
            bcj.ϵ__   .+= Δϵ
        else                    # plastic
            # Radial Return
            Δγ          = flow_rule / (2μ + 2(h + H) / 3)     # original
            n           = bcj.ξ__ ./ ξ_mag
            σ__prev     = bcj.σ__
            bcj.σ__    .= σₜᵣ__ - (2μ * Δγ) .* n
            bcj.α__    .= αₜᵣ__ + ( h * Δγ) .* n
            bcj.κ       = κₜᵣ   + (H * Δγ)  # original
            bcj.ϵₚ__  .+= (Δϵ - ((bcj.σ__ - σ__prev) ./ 2μ))
            bcj.ϵ__   .+= Δϵ
        end
        bcj.ϵ_dot_plastic__ .= (f * sinh(V \ (ξ_mag - bcj.κ - Y)) / ξ_mag) .* bcj.ξ__
        record!(history, i, bcj)
    end
    return nothing
end

function solve!(bcj::BCJMetalCurrentConfiguration{DK, <:AbstractFloat},
        history::BCJMetalConfigurationHistory)
    μ, Δϵ, Δt       = bcj.μ, bcj.Δϵ, bcj.Δt
    ϵ_dot_effective = bcj.ϵ_dot_effective
    V, Y, f, β      = bcj.V, bcj.Y, bcj.f, bcj.β
    h, r_d, r_s     = bcj.h, bcj.r_d, bcj.r_s # alpha
    H, R_d, R_s     = bcj.H, bcj.R_d, bcj.R_s # kappa
    sqrt23          = √(2 / 3)
    # timestep calculations
    for i ∈ range(2, bcj.N + 1)
        α_mag       = symmetricmagnitude(bcj.α__)
        # α_mag = sqrt( α_mag * 3./2.)       # match cho
        α_mag      *= sqrt23       # match vumat20
        # trial guesses: ISVs (from recovery) and stress
        recovery    = Δt * (r_d * ϵ_dot_effective + r_s) * α_mag    # recovery for alpha (kinematic hardening)
        Recovery    = Δt * (R_d * ϵ_dot_effective + R_s) * bcj.κ    # recovery for kappa (isotropic hardening)
        αₜᵣ__       = bcj.α__  .* (1 - recovery)
        κₜᵣ         = bcj.κ     * (1 - Recovery)

        ## trial stress guess
        σₜᵣ__       = bcj.σ__ + 2μ .* Δϵ           # trial stress
        bcj.ξ__    .= σₜᵣ__ - (2. / 3.) .* αₜᵣ__       # trial overstress original
        # ξ__          .= σₜᵣ__ - sqrt23 .* αₜᵣ__   # trial overstress FIT
        ξ_mag       = symmetricmagnitude(bcj.ξ__)



        # ----------------------------------- #
        ###   ---   YIELD CRITERION   ---   ###
        # ----------------------------------- #
        flow_rule = ξ_mag - sqrt23 * (κₜᵣ + β)         # same as vumat20
        # Crit = Xi_mag - (Katr + β) #changed to FIT
        if flow_rule <= 0.      # elastic
            # trial guesses are correct
            bcj.σ__    .= σₜᵣ__
            bcj.α__    .= αₜᵣ__
            bcj.κ       = κₜᵣ
            bcj.ϵ__   .+= Δϵ
        else                    # plastic
            # Radial Return
            Δγ          = flow_rule / (2μ + 2(h + H) / 3)     # original
            n           = bcj.ξ__ ./ ξ_mag
            σ__prev     = bcj.σ__
            bcj.σ__    .= σₜᵣ__ - (2μ * Δγ) .* n
            bcj.α__    .= αₜᵣ__ + ( h * Δγ) .* n
            bcj.κ       = κₜᵣ   + (H * sqrt23 * Δγ)  # original
            bcj.ϵₚ__  .+= (Δϵ - ((bcj.σ__ - σ__prev) ./ 2μ))
            bcj.ϵ__   .+= Δϵ
        end
        bcj.ϵ_dot_plastic__ .= (f * sinh(V \ (ξ_mag - bcj.κ - Y)) / ξ_mag) .* bcj.ξ__
        record!(history, i, bcj)
    end
    return nothing
end