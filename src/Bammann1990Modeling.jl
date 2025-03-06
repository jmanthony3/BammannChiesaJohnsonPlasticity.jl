# using PlasticityBase
using ComponentArrays
using ContinuumMechanicsBase
import ForwardDiff
using LossFunctions
using Optimization
using StructArrays
using Tensors

export Bammann1990Modeling, update, AbstractBCJMetalTest, predict, BCJMetalDataEntry, BCJMetalUniaxialTest, parameters, parameter_bounds, BCJProblem

# mutable struct Bammann1990Modeling{T, S<:SymmetricTensor{2, 3, T}} <: BCJMetal
struct Bammann1990Modeling <: BCJMetal
    θ               # ::T         # applied temperature
    ϵ_dot_effective # ::T         # strain rate (effective)
    ϵₙ              # ::T         # final strain
    μ               # ::T         # shear modulus at temperature, θ
    N               # ::Integer   # number of strain increments
    Δϵ              # ::Vector{T} # S         # total strain tensor step
    Δt              # ::T         # time step
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
    # # array declarations
    # ## OSVs
    # σ__             = zero(S) # zeros(T, 6)       # deviatoric stress
    # ϵₚ__            = zero(S) # zeros(T, 6)       # plastic strain
    # ϵ_dot_plastic__ = zero(S) # zeros(T, 6)       # plastic strain rate
    # ϵ__             = zero(S) # zeros(T, 6)       # total strain
    # ## ISVs
    # α__             = fill(1e-7, S) # fill(1e-7, 6) # alpha: kinematic hardening
    # κ               = 0.            # kappa: isotropic hardening
    # ## holding values
    Δϵ              = zero(S) # zeros(T, 6)       # strain increment
    # ξ__             = zero(S) # zeros(T, 6)       # overstress (S - 2/3*alpha)
    Δt              = (ϵₙ / N) / ϵ_dot

    # state evaluation - loading type
    ϵ_dot_effective = if loadtype ∈ (:tension, :compression)    # uniaxial tension/compression
        δϵ  = ϵₙ / N
        # Δϵ .= [δϵ, -0.499δϵ, -0.499δϵ, 0., 0., 0.]
        Δϵ  = S([δϵ, 0.0, 0.0, -0.499δϵ, 0.0, -0.499δϵ])
        Δt  = δϵ / ϵ_dot            # timestep
        ϵ_dot
    elseif loadtype == :torsion                                 # torsion
        # convert equivalent strain to true shear strain
        ϵₙ *= 0.5 * √(3.)
        # Δϵ .= [0., 0., 0., ϵₙ / N, 0., 0.]
        Δϵ  = S([0.0, ϵₙ / N, 0.0, 0.0, 0.0, 0.0])
        # equivalent strain rate to true shear strain rate
        Δt  = Δϵ[1, 2] / ϵ_dot      # timestep
        2ϵ_dot / √3.
    end
    return Bammann1990Modeling(θ, ϵ_dot_effective, ϵₙ, μ, N, Δϵ, Δt)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function update(model::Bammann1990Modeling, σ__, α__, κ, ϵ__, ϵₚ__, (;
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
    θ       = model.θ
    ϵ_dot_effective   = model.ϵ_dot_effective
    μ       = model.μ
    ϵₙ      = model.ϵₙ
    N       = model.N
    Δϵ      = model.Δϵ
    Δt      = model.Δt
    # loadtype= model.loadtype
    M       = N + 1
    T       = typeof(float(θ))
    S       = SymmetricTensor{2, 3, T}
    # # array declarations
    # ## OSVs
    # σ__             = zero(S) # zeros(T, 6)       # deviatoric stress
    # ϵₚ__            = zero(S) # zeros(T, 6)       # plastic strain
    # ϵ_dot_plastic__ = zero(S) # zeros(T, 6)       # plastic strain rate
    # ϵ__             = zero(S) # zeros(T, 6)       # total strain
    # ## ISVs
    # α__             = fill(1e-7, S) # fill(1e-7, 6) # alpha: kinematic hardening
    # κ               = 0.            # kappa: isotropic hardening
    # ## holding values
    # Δϵ              = zero(S) # zeros(T, 6)       # strain increment
    # ξ__             = zero(S) # zeros(T, 6)       # overstress (S - 2/3*alpha)
    # T = typeof(θ)
    # temperature dependent constants
    V   = C₁    * exp( -C₂ / θ )
    Y   = C₃    * exp(  C₄ / θ )
    f   = C₅    * exp( -C₆ / θ )
    β   = Y + (V * asinh( ϵ_dot_effective / f ))
    r_d = C₇    * exp( -C₈  / θ )
    h   = C₉    * exp(  C₁₀ * θ )
    r_s = C₁₁   * exp( -C₁₂ / θ )
    R_d = C₁₃   * exp( -C₁₄ / θ )
    H   = C₁₅   * exp(  C₁₆ * θ )
    R_s = C₁₇   * exp( -C₁₈ / θ )

    # # timestep calculations
    # for i ∈ range(2, M)
    # end
    # trial guesses
    α_mag       = norm(α__)
    # trial guesses: ISVs (from recovery) and stress
    recovery    = Δt * (r_d * ϵ_dot_effective + r_s) * α_mag    # recovery for alpha (kinematic hardening)
    Recovery    = Δt * (R_d * ϵ_dot_effective + R_s) * κ  # recovery for kappa (isotropic hardening)
    # αₜᵣ__       = state.α__  .* (1 - recovery)
    αₜᵣ__       = α__   * (1 - recovery)
    κₜᵣ         = κ     * (1 - Recovery)
    σₜᵣ__       = σ__ + (2μ * Δϵ)                         # trial stress
    ξ__   = @. σₜᵣ__ - αₜᵣ__                                             # trial overstress original
    # ξ__         = σₜᵣ__ - αₜᵣ__                                             # trial overstress original
    # ξ__          .= σₜᵣ__ - sqrt23 .* αₜᵣ__                                 # trial overstress FIT
    ξ_mag       = norm(ξ__)

    # yield criterion
    flow_rule = ξ_mag - κₜᵣ - β                                             # same as vumat20
    if flow_rule <= 0.      # elastic
        # trial guesses are correct
        σ__               = @. σₜᵣ__
        α__               = @. αₜᵣ__
        κ                 = κₜᵣ
        # state.ξ__               = ξ__
        ϵ__              += @. Δϵ
        # state.ϵ_dot_plastic__    .= 0.
    else                    # plastic
        # Radial Return
        Δγ                      = flow_rule / (2μ + 2(h + H) / 3)     # original
        # n̂                       = state.ξ__ ./ ξ_mag
        n̂                       = ξ__ / ξ_mag
        σ__prev                 = σ__
        # state.σ__               = @. σₜᵣ__ - (2state.μ * Δγ) .* n̂
        # state.α__               = @. αₜᵣ__ + ( h * Δγ) .* n̂
        σ__               = @. σₜᵣ__ - (2μ * Δγ) * n̂
        α__               = @. αₜᵣ__ + ( h * Δγ) * n̂
        # state.ξ__               = state.σ__ - state.α__
        κ                 = κₜᵣ   + (H * Δγ)  # original
        # state.ϵₚ__             += @. (state.Δϵ - ((state.σ__ - σ__prev) ./ 2state.μ))
        ϵₚ__             += @. (Δϵ - ((σ__ - σ__prev) / μ))
        ϵ__              += @. Δϵ
    end
    ϵ_dot_plastic__  = @. (f * sinh(V \ (ξ_mag - κ - Y)) / ξ_mag) * ξ__
    return σ__, α__, κ, ϵ__, ϵₚ__
end

# [20250301T0101] (JMA3)
    # - so far, the above is tested to march through the number of strain increments, N
    # - I suspect, that for the below to work, the above must just take the current material state
    #       and applied strain increment to return the deviatoric stress for the below to work.

abstract type AbstractBCJMetalTest{T, S} <: ContinuumMechanicsBase.AbstractMaterialTest end

## Predict overloads
function ContinuumMechanicsBase.predict(
            ψ   ::Bammann1990Modeling,#{T, S},
            test::AbstractBCJMetalTest{T, T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat, S<:SymmetricTensor{2, 3, T}}
    ϕ = ψ
    # ϕ = deepcopy(ψ)
    s = SymmetricTensor{2, 3, Float64}
    # if !isnothing(ui)
    #     for (name, value) in zip(keys(ui), ui)
    #         p[name] = value
    #     end
    # end
    # ϵ⃗ = []
    # σ⃗ = []
    # push!(ϵ⃗,ψ.ϵ__)
    # push!(σ⃗, ψ.σ__)
    # ϵ⃗ = zeros(T, (6, ψ.N + 1)) # zeros(S, ψ.N + 1) # ψ.ϵ__ .+ [ψ.Δϵ * i for i ∈ range(0, test.N)]
    # σ⃗ = zeros(T, (6, ψ.N + 1)) # zeros(S, ψ.N + 1)
    # ϵ⃗[:, 1], σ⃗[:, 1] = ϕ.ϵ__, ϕ.σ__
    θ       = ϕ.θ
    ϵ_dot_effective   = ϕ.ϵ_dot_effective
    μ       = ϕ.μ
    ϵₙ      = ϕ.ϵₙ
    N       = ϕ.N
    Δϵ      = ϕ.Δϵ
    Δt      = ϕ.Δt
    # loadtype= ϕ.loadtype
    M       = N + 1
    # T       = typeof(float(θ))
    # S       = SymmetricTensor{2, 3, T}
    # array declarations
    ## OSVs
    σ__             = zero(s) # zeros(T, 6)       # deviatoric stress
    ϵₚ__            = zero(s) # zeros(T, 6)       # plastic strain
    ϵ_dot_plastic__ = zero(s) # zeros(T, 6)       # plastic strain rate
    ϵ__             = zero(s) # zeros(T, 6)       # total strain
    ## ISVs
    α__             = fill(1e-7, s) # fill(1e-7, 6) # alpha: kinematic hardening
    κ               = 0.            # kappa: isotropic hardening
    ## holding values
    Δϵ              = zero(s) # zeros(T, 6)       # strain increment
    ξ__             = zero(s) # zeros(T, 6)       # overstress (S - 2/3*alpha)
    ϵ⃗ = zeros(s, M) # ψ.ϵ__ .+ [ψ.Δϵ * i for i ∈ range(0, test.N)]
    σ⃗ = zeros(s, M)
    ϵ⃗[1], σ⃗[1] = s(ϵ__), s(σ__)
    # @show α__, κ
    for i ∈ range(2, M)
        σ__, α__, κ, ϵ__, ϵₚ__ = update(ϕ, σ__, α__, κ, ϵ__, ϵₚ__, p)
        # @show α__, κ
        # push!(ϵ⃗, ϕ.ϵ__)
        # push!(σ⃗, ϕ.σ__)
        # ϵ⃗[:, i], σ⃗[:, i] = ϕ.ϵ__, ϕ.σ__
        # @show s(ϕ.ϵ__)
        # @show s(ϕ.σ__)
        ϵ⃗[i], σ⃗[i] = s(ϵ__), s(σ__)
    end
    return (data=(λ=ϵ⃗, s=σ⃗),)
    # return (data=(λ=hcat(ϵ⃗...), s=hcat(σ⃗...)),)
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
    ψ   ::Bammann1990Modeling,#{T, S},
    test::BCJMetalUniaxialTest{T, T},
    u0;
    ad_type,
    ui,
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
        ψ, test, qs, loss, ad_type, kwargs = p
        function g(qs)
            for (name, value) in zip(keys(qs), qs)
                if !isnan(value)
                    # @show value
                    ps[name] = value
                end
            end
            return ComponentVector(ps)
        end
        ϕ = deepcopy(ψ)
        # @show g(qs)
        pred = predict(ϕ, test, g(qs); ad_type, kwargs...)
        res = map(i -> loss.(symmetricvonMises(i[1]), only(i[2])), zip(pred.data.s[2:end], test.data.s)) |> mean
        @show res
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
    # ϕ = ψ
    ϕ = deepcopy(ψ)
    p = (ϕ, test, ui, loss, ad_type, kwargs)
    return OptimizationProblem(func, u0, p; lb, ub, int, lcons, ucons, sense)
end