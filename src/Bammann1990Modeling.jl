# using PlasticityBase
using ComponentArrays
using ContinuumMechanicsBase
import ForwardDiff
using LossFunctions
using Optimization
using StructArrays
using Tensors

export Bammann1990Modeling, update!, AbstractBCJMetalTest, predict, BCJMetalDataEntry, BCJMetalUniaxialTest, parameters, parameter_bounds, BCJProblem

mutable struct Bammann1990Modeling{T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}} <: BCJMetal
    N               ::Integer   # number of strain increments
    θ               ::T         # applied temperature
    μ               ::T         # shear modulus at temperature, θ
    σ__             ::S         # deviatoric stress tensor
    ϵₚ__            ::S         # plastic strain tensor
    ϵ_dot_plastic__ ::S         # plastic strain rate
    ϵ__             ::S         # total strain tensor
    ϵ_dot_effective ::T         # strain rate (effective)
    Δϵ              ::S         # total strain tensor step
    Δt              ::T         # time step
    α__             ::S         # kinematic hardening tensor
    κ               ::T         # isotropic hardening scalar
    ξ__             ::S         # overstress tensor (S - 2/3*alpha)
end

function Bammann1990Modeling(bcj::BCJMetalStrainControl, μ::AbstractFloat)
    θ       = bcj.θ
    ϵ_dot   = bcj.ϵ_dot
    ϵₙ      = bcj.ϵₙ
    N       = bcj.N
    loadtype= bcj.loadtype
    M       = N + 1
    T       = typeof(float(θ))
    S       = SymmetricTensor{2, 3, T}
    # array declarations
    ## OSVs
    σ__             = zero(S)       # deviatoric stress
    ϵₚ__            = zero(S)       # plastic strain
    ϵ_dot_plastic__ = zero(S)       # plastic strain rate
    ϵ__             = zero(S)       # total strain
    ## ISVs
    α__             = fill(1e-7, S) # alpha: kinematic hardening
    κ               = 0.            # kappa: isotropic hardening
    ## holding values
    Δϵ              = zero(S)       # strain increment
    ξ__             = zero(S)       # overstress (S - 2/3*alpha)

    # state evaluation - loading type
    ϵ_dot_effective = if loadtype ∈ (:tension, :compression)    # uniaxial tension/compression
        δϵ  = ϵₙ / N
        Δϵ  = S([δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ])
        Δt  = δϵ / ϵ_dot            # timestep
        ϵ_dot
    elseif loadtype == :torsion                                 # torsion
        # convert equivalent strain to true shear strain
        ϵₙ *= 0.5 * √(3.)
        Δϵ  = S([0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0])
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[1, 2] / ϵ_dot      # timestep
        2ϵ_dot / √3.
    end
    return Bammann1990Modeling{T, S}(N, θ, μ,
        σ__, ϵₚ__, ϵ_dot_plastic__, ϵ__, ϵ_dot_effective, Δϵ, Δt, α__, κ, ξ__)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function update!(state::Bammann1990Modeling, (;
            C₁,     C₂,     # V
            C₃,     C₄,     # Y
            C₅,     C₆,     # f
            C₇,     C₈,     # r_d
            C₉,     C₁₀,    # h
            C₁₁,    C₁₂,    # r_s
            C₁₃,    C₁₄,    # R_d
            C₁₅,    C₁₆,    # H
            C₁₇,    C₁₈     # R_s
        ))
    T = typeof(state.θ)
    # temperature dependent constants
    V   = C₁    * exp( -C₂ / state.θ )
    Y   = C₃    * exp(  C₄ / state.θ )
    f   = C₅    * exp( -C₆ / state.θ )
    β   = Y + (V * asinh( state.ϵ_dot_effective / f ))
    r_d = C₇    * exp( -C₈  / state.θ )
    h   = C₉    * exp(  C₁₀ * state.θ )
    r_s = C₁₁   * exp( -C₁₂ / state.θ )
    R_d = C₁₃   * exp( -C₁₄ / state.θ )
    H   = C₁₅   * exp(  C₁₆ * state.θ )
    R_s = C₁₇   * exp( -C₁₈ / state.θ )

    # # timestep calculations
    # for i ∈ range(2, bcj.N + 1)
    # end
    # trial guesses
    α_mag       = norm(state.α__)
    # trial guesses: ISVs (from recovery) and stress
    recovery    = state.Δt * (r_d * state.ϵ_dot_effective + r_s) * α_mag    # recovery for alpha (kinematic hardening)
    Recovery    = state.Δt * (R_d * state.ϵ_dot_effective + R_s) * state.κ  # recovery for kappa (isotropic hardening)
    αₜᵣ__       = state.α__   * (1 - recovery)
    κₜᵣ         = state.κ     * (1 - Recovery)
    σₜᵣ__       = state.σ__ + (2state.μ * state.Δϵ)                         # trial stress
    ξ__         = σₜᵣ__ - αₜᵣ__                                             # trial overstress original
    # ξ__          .= σₜᵣ__ - sqrt23 .* αₜᵣ__                                 # trial overstress FIT
    ξ_mag       = norm(state.ξ__)

    # yield criterion
    flow_rule = ξ_mag - κₜᵣ - β                                             # same as vumat20
    if flow_rule <= 0.      # elastic
        # trial guesses are correct
        state.σ__               = σₜᵣ__
        state.α__               = αₜᵣ__
        state.κ                 = κₜᵣ
        state.ξ__               = ξ__
        state.ϵ__              += state.Δϵ
        # state.ϵ_dot_plastic__    .= 0.
    else                    # plastic
        # Radial Return
        Δγ                      = flow_rule / (2state.μ + 2(h + H) / 3)     # original
        n̂                       = state.ξ__ / ξ_mag
        σ__prev                 = state.σ__
        state.σ__               = σₜᵣ__ - (2state.μ * Δγ) * n̂
        state.α__               = αₜᵣ__ + ( h * Δγ) * n̂
        state.ξ__               = state.σ__ - state.α__
        state.κ                 = κₜᵣ   + (H * Δγ)  # original
        state.ϵₚ__             += (state.Δϵ - ((state.σ__ - σ__prev) / 2state.μ))
        state.ϵ__              += state.Δϵ
    end
    state.ϵ_dot_plastic__   = (f * sinh(V \ (ξ_mag - state.κ - Y)) / ξ_mag) * state.ξ__
    return nothing
end

# [20250301T0101] (JMA3)
    # - so far, the above is tested to march through the number of strain increments, N
    # - I suspect, that for the below to work, the above must just take the current material state
    #       and applied strain increment to return the deviatoric stress for the below to work.

abstract type AbstractBCJMetalTest{T, S} <: ContinuumMechanicsBase.AbstractMaterialTest end

## Predict overloads
function ContinuumMechanicsBase.predict(
            ψ   ::Bammann1990Modeling{T, S},
            test::AbstractBCJMetalTest{T, T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}}
    ϵ⃗ = zeros(S, ψ.N + 1) # ψ.ϵ__ .+ [ψ.Δϵ * i for i ∈ range(0, test.N)]
    σ⃗ = zeros(S, ψ.N + 1)
    ϵ⃗[1], σ⃗[1] = ψ.ϵ__, ψ.σ__
    for i ∈ range(2, ψ.N + 1)
        update!(ψ, p)
        ϵ⃗[i], σ⃗[i] = ψ.ϵ__, ψ.σ__
    end
    return (data=(λ=ϵ⃗, s=σ⃗),)
end

struct BCJMetalDataEntry{T, S}
    λ::Vector{T}
    s::Vector{S}
end

struct BCJMetalUniaxialTest{T, S} <: AbstractBCJMetalTest{T, S}
    data::StructVector
    name::String
    """
    $(SIGNATURES)

    Creates an object storing results from a uniaxial test of a hyperelatic  material.

    # Arguments:
    - `λ₁`: Vector of uniaxial stretches
    - `s₁`: Vector of experiemntal stresses (optional)
    - `name`: string for the name of the test
    - `incompressible`: `true` if the material can be assumed to be incompressible.
    """
    function BCJMetalUniaxialTest(λ₁, s₁; name, incompressible = true)
        @assert length(λ₁) == length(s₁) "Inputs must be the same length"
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁))
        data = StructArray{BCJMetalDataEntry}((λ, s))
        new{eltype(eltype(λ)),eltype(eltype(s))}(data, name)
    end
    function BCJMetalUniaxialTest(λ₁; name, incompressible = true)
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(Vector{eltype(λ₁)}(undef, length(λ₁))))
        data = StructArray{BCJMetalDataEntry}((λ, s))
        new{eltype(eltype(λ)),eltype(eltype(s))}(data, name)
    end
end

parameters(::Bammann1990Modeling) = (
    :C₁,    :C₂,     # V
    :C₃,    :C₄,     # Y
    :C₅,    :C₆,     # f
    :C₇,    :C₈,     # r_d
    :C₉,    :C₁₀,    # h
    :C₁₁,   :C₁₂,    # r_s
    :C₁₃,   :C₁₄,    # R_d
    :C₁₅,   :C₁₆,    # H
    :C₁₇,   :C₁₈     # R_s
)

function parameter_bounds(::ContinuumMechanicsBase.AbstractMaterialModel, ::Any)
    lb = nothing
    ub = nothing
    return (lb = lb, ub = ub)
end

function parameter_bounds(
            ψ       ::ContinuumMechanicsBase.AbstractMaterialModel,
            tests   ::Vector{Any},
        )
    bounds = map(Base.Fix1(parameter_bounds, ψ), tests)
    lbs = getfield.(bounds, :lb)
    ubs = getfield.(bounds, :ub)
    if !(eltype(lbs) <: Nothing)
        lb_ps = fieldnames(eltype(lbs))
        lb = map(p -> p .=> maximum(getfield.(lbs, p)), lb_ps) |> NamedTuple
    else
        lb = nothing
    end
    if !(eltype(ubs) <: Nothing)
        ub_ps = fieldnames(eltype(ubs))
        ub = map(p -> p .=> minimum(getfield.(ubs, p)), ub_ps) |> NamedTuple
    else
        ub = nothing
    end
    return (lb = lb, ub = ub)
end

function BCJProblem(
    ψ   ::Bammann1990Modeling{T, S},
    test::BCJMetalUniaxialTest{T, T},
    u0;
    ad_type,
    loss    = L2DistLoss(),
    lb      = parameter_bounds(ψ, test).lb,
    ub      = parameter_bounds(ψ, test).ub,
    int     = nothing,
    lcons   = nothing,
    ucons   = nothing,
    sense   = nothing,
    kwargs...,
) where {T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}}
    function f(ps, p)
        ψ, test, loss, ad_type, kwargs = p
        pred = predict(ψ, test, ps; ad_type, kwargs...)
        res = map(i -> loss.(i[1], i[2]), zip(pred.data.s, test.data.s)) |> mean
        return res
    end

    u0 = ComponentVector(u0)
    if !isnothing(lb) && !isnothing(ub)
        lb = ComponentVector(lb)
        ub = ComponentVector(ub)
    elseif !isnothing(lb)
        lb = ComponentVector(lb)
        ub = u0 .* Inf
    elseif !isnothing(ub)
        ub = ComponentVector(ub)
        lb = u0 .* -Inf
    else
        ub = u0 .* Inf
        lb = u0 .* -Inf
    end

    model_ps = parameters(ψ)
    for p in model_ps
        if !isnothing(lb)
            if (u0[p] < lb[p])
                @error "Parameter $p = $(u0[p]) is less than lower bound of $(lb[p])"
                return nothing
            end
        end
        if !isnothing(ub)
            if (u0[p] > ub[p])
                @error "Parameter $p = $(u0[p]) is greater than upper bound of $(ub[p])"
                return nothing
            end
        end
    end

    func = OptimizationFunction(f, ad_type)
    # Check for Bounds
    p = (ψ, test, loss, ad_type, kwargs)
    return OptimizationProblem(func, u0, p; lb, ub, int, lcons, ucons, sense)
end