ContinuumMechanicsBase.parameters(::Bammann1990Modeling) = (
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

function ContinuumMechanicsBase.parameter_bounds(::Bammann1990Modeling, ::T) where {T<:ContinuumMechanicsBase.AbstractMaterialTest}
    lb = (
            C₁  = 0.0,  C₂  = 0.0,  # V
            C₃  = 0.0,  C₄  = 0.0,  # Y
            C₅  = 0.0,  C₆  = 0.0,  # f
            C₇  = 0.0,  C₈  = 0.0,  # r_d
            C₉  = 0.0,  C₁₀ = 0.0,  # h
            C₁₁ = 0.0,  C₁₂ = 0.0,  # r_s
            C₁₃ = 0.0,  C₁₄ = 0.0,  # R_d
            C₁₅ = 0.0,  C₁₆ = 0.0,  # H
            C₁₇ = 0.0,  C₁₈ = 0.0,  # R_s
        )
    ub = nothing
    return (lb = lb, ub = ub)
end

function ContinuumMechanicsBase.MaterialOptimizationProblem(
    ψ   ::Bammann1990Modeling{T},  # , S},
    test::BCJMetalUniaxialTest{T},
    u₀,
    ps,
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
        pred = ContinuumMechanicsBase.predict(ψ, test, g(ps, qs); ad_type, kwargs...)
        resϵ = [first(x) for x in eachcol(pred.data.ϵ)]
        testϵ = [first(x) for x in test.data.ϵ]
        # resϵ = [x[1, 1] for x in pred.data.ϵ]
        # testϵ = [x[1, 1] for x in test.data.ϵ]
        s = collect([[x...] for x in eachcol(pred.data.σ)[[findlast(x .>= resϵ) for x in testϵ]]])
        # s = collect([[x...] for x in pred.data.σ[[findlast(x .>= resϵ) for x in testϵ]]])
        res = map(i -> loss.(symmetricvonMises(i[1]), only(i[2])), zip(s, test.data.σ)) |> mean
        @show res
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

    model_ps = ContinuumMechanicsBase.parameters(ψ)
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