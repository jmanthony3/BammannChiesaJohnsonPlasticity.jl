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
    get_data(d, f; columnate=false, concatenate=false) = (z = [[first(x) for x in (columnate ? eachcol(y.data[f]) : y.data[f])] for y in d]; concatenate ? vcat(z...) : z)
    get_idx(d1, d2) = [[findlast(first(x) .>= get_data([d2[i]], :ϵ; columnate=true)...) for x in z] for (i, z) in enumerate([[first(x) for x in y.data.ϵ] for y in values(d1)])]
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
        preds = [ContinuumMechanicsBase.predict(ψ, test, g(ps, qs); ad_type, kwargs...) for (ψ, test) in zip(ψs, tests)]
        # @show preds
        # resϵ = [first(x) for x in eachcol(pred.data.ϵ)]
        # testϵ = [first(x) for x in test.data.ϵ]
        # ϵ = get_data(tests, ϵ)
        # ϵ̂ = get_data(preds, ϵ; col=true)
        # # resϵ = [x[1, 1] for x in pred.data.ϵ]
        # # testϵ = [x[1, 1] for x in test.data.ϵ]
        # s = collect([[x...] for x in eachcol(pred.data.σ)[[findlast(x .>= resϵ) for x in testϵ]]])
        # # s = collect([[x...] for x in pred.data.σ[[findlast(x .>= resϵ) for x in testϵ]]])
        s = vcat([x[y] for (x, y) in zip(get_data(preds, :ϵ; columnate=true), get_idx(tests, preds))]...)
        ŝ = vcat([[first(x) for x in y.data.ϵ] for y in values(tests)]...)
        # @show length(s) == length(ŝ)
        res = map(i -> loss.(only(i[1]), only(i[2])), zip(s, ŝ)) |> mean
        # @show res # uncomment for testing
        return res
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