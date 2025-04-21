module OptimizationBCJPlasticityExt

using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays
using Optimization, Interpolations, LossFunctions

export parameter_bounds, MaterialOptimizationProblem

"Set lower bounds to zero for physical admissibility."
function ContinuumMechanicsBase.parameter_bounds(ψ::AbstractBCJModel, test::AbstractBCJTest)
    lb = NamedTuple(Symbol.(parameters(ψ)) .=> 0.0)
    ub = nothing
    return (lb = lb, ub = ub)
end

"Dispatch for BCJ-specific types and functions."
function ContinuumMechanicsBase.MaterialOptimizationProblem(
    ψ   ::AbstractBCJModel,  # , S},
    test::AbstractBCJTest,
    u₀,
    model_ps,
    ad_type,
    loss;
    ui,
    lb      = parameter_bounds(ψ, test).lb,
    ub      = parameter_bounds(ψ, test).ub,
    int     = nothing,
    lcons   = nothing,
    ucons   = nothing,
    sense   = nothing,
    kwargs...,
) # where {T<:AbstractFloat} #, S<:SymmetricTensor{2, 3, T}}
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
        prediction = ContinuumMechanicsBase.predict(ψ, test, g(ps, qs); ad_type, kwargs...)
        ϵ = [first(x) for x in test.data.ϵ]
        σ = [first(x) for x in test.data.σ]
        ϵ̂ = [first(x) for x in eachcol(prediction.data.ϵ)]
        σ̂ = collect(eachcol(prediction.data.σ))
        ŝ = linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵ̂)
        err = map(i -> loss.(i[1], vonMises(i[2])), zip(ŝ, σ̂)) |> mean
        return err
    end

    u₀ = ComponentVector(u₀)
    # pb = ContinuumMechanicsBase.parameter_bounds(ψ, test)
    # lb, ub = pb.lb, pb.ub
    if !isnothing(lb) && !isnothing(ub)
        lb = ComponentVector(lb)
        ub = ComponentVector(ub)
    elseif !isnothing(lb)
        lb = ComponentVector(lb)
        ub = u₀ .* Inf
    elseif !isnothing(ub)
        ub = ComponentVector(ub)
        lb = u₀ .* -Inf
    else
        ub = u₀ .* Inf
        lb = u₀ .* -Inf
    end

    # model_ps = ContinuumMechanicsBase.parameters(ψ)
    for p in model_ps
        if !isnothing(lb)
            if (u₀[p] < lb[p])
                @error "Parameter $p = $(u₀[p]) is less than lower bound of $(lb[p])"
                return nothing
            end
        end
        if !isnothing(ub)
            if (u₀[p] > ub[p])
                @error "Parameter $p = $(u₀[p]) is greater than upper bound of $(ub[p])"
                return nothing
            end
        end
    end

    func = OptimizationFunction(f, ad_type)
    # Check for Bounds
    p = (ψ, test, ui, loss, ad_type, kwargs)
    return OptimizationProblem(func, u₀, p; lb, ub, int, lcons, ucons, sense)
end

function ContinuumMechanicsBase.MaterialOptimizationProblem(
    ψs   ::Vector{<:AbstractBCJModel},  # , S},
    tests::Vector{<:AbstractBCJTest},
    u₀,
    model_ps,
    ad_type,
    loss;
    ui,
    lb      = parameter_bounds(first(ψs), first(tests)).lb,
    ub      = parameter_bounds(first(ψs), first(tests)).ub,
    int     = nothing,
    lcons   = nothing,
    ucons   = nothing,
    sense   = nothing,
    kwargs...,
) # where {T<:AbstractFloat} #, S<:SymmetricTensor{2, 3, T}}
    function f(ps, p)
        ψs, tests, qs, loss, ad_type, kwargs = p
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
        errors = Vector{typeof(first(ψs).θ)}(undef, length(ψs))
        for (i, (ψ, test)) in enumerate(zip(ψs, tests))
            prediction = ContinuumMechanicsBase.predict(ψ, test, g(ps, qs); ad_type, kwargs...)
            ϵ = [first(x) for x in test.data.ϵ]
            σ = [first(x) for x in test.data.σ]
            ϵ̂ = [first(x) for x in eachcol(prediction.data.ϵ)]
            σ̂ = collect(eachcol(prediction.data.σ))
            ŝ = linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵ̂)
            errors[i] = map(i -> loss.(i[1], vonMises(i[2])), zip(ŝ, σ̂)) |> mean
        end
        return mean(errors)
    end

    u₀ = ComponentVector(u₀)
    # pb = ContinuumMechanicsBase.parameter_bounds(ψ, test)
    # lb, ub = pb.lb, pb.ub
    if !isnothing(lb) && !isnothing(ub)
        lb = ComponentVector(lb)
        ub = ComponentVector(ub)
    elseif !isnothing(lb)
        lb = ComponentVector(lb)
        ub = u₀ .* Inf
    elseif !isnothing(ub)
        ub = ComponentVector(ub)
        lb = u₀ .* -Inf
    else
        ub = u₀ .* Inf
        lb = u₀ .* -Inf
    end

    # model_ps = ContinuumMechanicsBase.parameters(ψ)
    for p in model_ps
        if !isnothing(lb)
            if (u₀[p] < lb[p])
                @error "Parameter $p = $(u₀[p]) is less than lower bound of $(lb[p])"
                return nothing
            end
        end
        if !isnothing(ub)
            if (u₀[p] > ub[p])
                @error "Parameter $p = $(u₀[p]) is greater than upper bound of $(ub[p])"
                return nothing
            end
        end
    end

    func = OptimizationFunction(f, ad_type)
    # Check for Bounds
    p = (ψs, tests, ui, loss, ad_type, kwargs)
    return OptimizationProblem(func, u₀, p; lb, ub, int, lcons, ucons, sense)
end

include("Metals.jl")

end # end of module
