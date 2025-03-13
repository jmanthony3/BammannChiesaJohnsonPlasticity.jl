using ContinuumMechanicsBase
using ComponentArrays, StructArrays
using Optimization, OptimizationOptimJL, LossFunctions
using DocStringExtensions

abstract type AbstractJCModel         <: ContinuumMechanicsBase.AbstractMaterialModel end
abstract type AbstractJCTest{T}       <: ContinuumMechanicsBase.AbstractMaterialTest end

struct JCStrainControl{T1<:Integer, T2<:AbstractFloat} <: ContinuumMechanicsBase.AbstractMaterialTest
    θ   ::T2    # applied temperature
    ϵ̇   ::T2    # applied strain rate
    ϵₙ  ::T2    # final strain
    N   ::T1    # number of strain increments
end

struct JCDataEntry{T, S}
    ϵ::Vector{T}
    σ::Vector{S}
end

struct JCUniaxialTest{T, S} <: AbstractJCTest{T}
    data::StructVector
    name::String
    """
    $(SIGNATURES)

    Creates an object storing results from a uniaxial test of a hyperelatic  material.

    # Arguments:
    - `ϵ₁`: Vector of experimental, uniaxial strains
    - `σ₁`: Vector of experimental, uniaxial stresses (optional)
    - `name`: string for the name of the test
    - `incompressible`: `true` if the material can be assumed to be incompressible.
    """
    function JCUniaxialTest(ϵ₁, σ₁; name, incompressible = true)
        @assert length(ϵ₁) == length(σ₁) "Inputs must be the same length"
        if incompressible
            # ϵ₂ = ϵ₃ = @. sqrt(1 / ϵ₁)
            ϵ₂ = ϵ₃ = @. -0.499ϵ₁
        else
            ϵ₂ = ϵ₃ = Vector{eltype(ϵ₁)}(undef, length(ϵ₁))
        end
        ϵ = collect.(zip(ϵ₁, ϵ₂, ϵ₃))
        σ = collect.(zip(σ₁))
        data = StructArray{JCDataEntry}((ϵ, σ))
        new{eltype(eltype(ϵ)), eltype(eltype(σ))}(data, name)
    end
    function JCUniaxialTest(ϵ₁; name, incompressible = true)
        if incompressible
            # ϵ₂ = ϵ₃ = @. sqrt(1 / ϵ₁)
            ϵ₂ = ϵ₃ = @. -0.499ϵ₁
        else
            ϵ₂ = ϵ₃ = Vector{eltype(ϵ₁)}(undef, length(ϵ₁))
        end
        ϵ = collect.(zip(ϵ₁, ϵ₂, ϵ₃))
        σ = collect.(zip(Vector{eltype(ϵ₁)}(undef, length(ϵ₁))))
        data = StructArray{JCDataEntry}((ϵ, σ))
        new{eltype(eltype(ϵ)), eltype(eltype(σ))}(data, name)
    end
end

struct JC{T<:AbstractFloat} <: AbstractJCModel
    θ   ::T         # applied temperature
    ϵ̇   ::T         # strain rate (effective)
    ϵₙ  ::T         # final strain
    N   ::Integer   # number of strain increments
    Tr  ::T         # reference temperature
    Tm  ::T         # melting temperature
    er0 ::T         # reference strain rate
    ϵ⁺  ::T
    θ⁺  ::T
    Δϵ  ::T         # strain increment
end

function JC(jc::JCStrainControl, Tr::T, Tm::T, er0::T) where {T<:AbstractFloat}
    ϵ⁺  = jc.ϵ̇ / er0
    θ⁺  = ( jc.θ - Tr ) / ( Tm - Tr )
    Δϵ  = jc.ϵₙ / jc.N # strain increment
    return JC{T}(jc.θ, jc.ϵ̇, jc.ϵₙ, jc.N, Tr, Tm, er0, ϵ⁺, θ⁺, Δϵ)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function update(model::JC, ϵ, (;
            A, B, n, C, m
        ))
    return (#=[=#   A   +   (  B * ( ϵ ^ n )  )     #=]=#) * (#=[=#
        1.0 + ( C * log(model.ϵ⁺) )                 #=]=#) * (#=[=#
        1.0 - ( model.θ⁺ ^ m )                      #=]=#)
end

function ContinuumMechanicsBase.predict(
            ψ   ::JC{T}, # , S},
            test::AbstractJCTest{T},
            p;
            kwargs...,
        ) where {T<:AbstractFloat} # , S<:SymmetricTensor{2, 3}}
    M = ψ.N + 1
    σ     = 0.0 # deviatoric stress
    ϵₚ    = 0.0 # plastic strain
    ϵ⃗ = []
    σ⃗ = []
    push!(ϵ⃗, ϵₚ)
    push!(σ⃗, σ)
    for i ∈ range(2, M)
        ϵₚ += ψ.Δϵ
        σ = update(ψ, ϵₚ, p)
        push!(ϵ⃗, ϵₚ)
        push!(σ⃗, σ)
    end
    return (data=(ϵ=hcat(ϵ⃗...), σ=hcat(σ⃗...)),)
end


function parameters(::M) where {M<:ContinuumMechanicsBase.AbstractMaterialModel} end

function parameter_bounds(::M, ::Any) where {M<:ContinuumMechanicsBase.AbstractMaterialModel}
    lb = nothing
    ub = nothing
    return (lb = lb, ub = ub)
end

function parameter_bounds(
            ψ       ::M,
            tests   ::Vector{Any},
        ) where {M<:ContinuumMechanicsBase.AbstractMaterialModel}
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

parameters(::JC) = (:A, :B, :n, :C, :m)

"""
$(SIGNATURES)

Creates an `OptimizationProblem` for use in [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) to find the optimal parameters.

# Arguments:
- `ψ`: material model to use
- `test` or `tests`: A single or vector of hyperelastics tests to use when fitting the parameters
- `u₀`: Initial guess for parameters
- `ps`: Any additional parameters for calling predict
- `adb`: Select differentiation type from [`ADTypes.jl`](https://github.com/SciML/ADTypes.jl). The type is automatically applied to the type of AD applied to the Optimization Problem also.
- `loss`: Loss function from [`LossFunctions.jl`](https://github.com/JuliaML/LossFunctions.jl)
"""
function JCPlasticityProblem end

function JCPlasticityProblem(
    ψ   ::JC{T}, # , S},
    test::JCUniaxialTest{T},
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
) where {T<:AbstractFloat} #, S<:SymmetricTensor{2, 3, T}}
    function f(ps, p)
        ψ, test, qs, loss, ad_type, kwargs = p
        function g(ps, qs)
            if !isnothing(qs) && any(!isnan, qs)
                for (name, value) in zip(keys(qs), qs)
                    if !isnan(value)
                        # @show value
                        ps[name] = value
                    end
                end
            end
            return ComponentVector(ps)
        end
        pred = predict(ψ, test, g(ps, qs); ad_type, kwargs...)
        resϵ = [first(x) for x in eachcol(pred.data.ϵ)]
        testϵ = [first(x) for x in test.data.ϵ]
        s = collect([[x...] for x in eachcol(pred.data.σ)[[findlast(x .>= resϵ) for x in testϵ]]])
        res = map(i -> loss.(only(i[1]), only(i[2])), zip(s, test.data.σ)) |> mean
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
    p = (ψ, test, ui, loss, ad_type, kwargs)
    return OptimizationProblem(func, u0, p; lb, ub, int, lcons, ucons, sense)
end