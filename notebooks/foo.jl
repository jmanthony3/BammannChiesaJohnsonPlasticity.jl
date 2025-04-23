using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, DataInterpolations, LossFunctions
using Plots, Printf

using Distributed, BenchmarkTools

include("Cho2019UnifiedStaticDynamic-functions.jl")

function MMO(
    ψ   ::Cho2019Unified,  # , S},
    test::BCJMetalUniaxialTest,
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
        # ŝ = linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵ̂)
        ŝ = CubicSpline(σ, ϵ; extrapolation=ExtrapolationType.Linear).(ϵ̂)
        err = map(i -> loss.(i[1], vonMises(i[2])), zip(ŝ, σ̂)) |> mean
        # ϵp = [first(x) for x in eachcol(prediction.data.ϵ)]
        # σp = collect(eachcol(prediction.data.σ))
        # # sp = linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵp)
        # sp = CubicSpline(σ, ϵ; extrapolation=ExtrapolationType.Linear).(ϵp)
        # err = map(i -> loss.(i[1], vonMises(i[2])), zip(sp, σp)) |> mean
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

function MMO(
    ψs   ::Vector{<:Cho2019Unified},  # , S},
    tests::Vector{<:BCJMetalUniaxialTest},
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
        # errors = Vector{typeof(first(ψs).θ)}(undef, length(ψs))
        errors = []
        for (i, (ψ, test)) in enumerate(zip(ψs, tests))
            prediction = ContinuumMechanicsBase.predict(ψ, test, g(ps, qs); ad_type, kwargs...)
            ϵ = [first(x) for x in test.data.ϵ]
            σ = [first(x) for x in test.data.σ]
            ϵ̂ = [first(x) for x in eachcol(prediction.data.ϵ)]
            σ̂ = collect(eachcol(prediction.data.σ))
            # ŝ = linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵ̂)
            ŝ = CubicSpline(σ, ϵ; extrapolation=ExtrapolationType.Linear).(ϵ̂)
            # errors[i] = map(i -> loss.(i[1], vonMises(i[2])), zip(ŝ, σ̂)) |> mean
            push!(errors, map(i -> loss.(i[1], vonMises(i[2])), zip(ŝ, σ̂)) |> mean)
            # ϵp = [first(x) for x in eachcol(prediction.data.ϵ)]
            # σp = collect(eachcol(prediction.data.σ))
            # # sp = linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵp)
            # sp = CubicSpline(σ, ϵ; extrapolation=ExtrapolationType.Linear).(ϵp)
            # push!(errors, map(i -> loss.(i[1], vonMises(i[2])), zip(sp, σp)) |> mean)
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

df_Fig4a = CSV.read("Cho2019UnifiedStaticDynamic-Fig4a.csv", DataFrame;
    header=true, delim=',', skipto=3, types=Float64)

n   = 2.0
ω₀  = 3.6e4
R   = 8.31446261815324 # universal gas constant
E⁺  = 82.0e3
z   = 0.65
d₀  = 10.0 # μm (Ghauri et al., 1990)
η₀  = 0.0
Kic = 1000.0
𝒹   = 0.0
𝒻   = 0.001
R₀  = 0.0

ϵ̇ = 4e-4

tests = Dict()
domains = Dict()
models = Dict()
for (i, θ) in enumerate((298, 407, 475, 509, 542, 559, 576, 610, 678, 814))
    θ_str = match(r"(.*)K(.*)", names(df_Fig4a)[4(i - 1) + 1])[1]
    θ_flt = parse(Float64, θ_str)
    x = filter(!ismissing, df_Fig4a[!, 4(i - 1) + 1])
    idx_sort = sortperm(x)
    x = x[idx_sort]
    y = filter(!ismissing, df_Fig4a[!, 4(i - 1) + 2])[idx_sort] .* 1e6
    @show (4(i - 1) + 1, 4(i - 1) + 2), θ_str, ϵ̇, last(x), 4length(x)
    tests[θ_str] = BCJMetalUniaxialTest(x, y, name="$(θ_flt)K")
    domains[θ_str] = BCJMetalStrainControl(θ_flt, ϵ̇, last(x), 4length(x), :tension)
    models[θ_str] = Cho2019Unified(domains[θ_str], n, ω₀, E⁺, E⁺, R, d₀, z, Kic, 𝒹, 𝒻, η₀, R₀)
end

tests = sort(tests; rev=false)
domains = sort(domains; rev=false)
models = sort(models; rev=false)

p = ComponentVector(
    C₁ = 5.637,
    C₂ = 112.6,
    C₃ = 8.378,
    C₄ = 324.9,
    C₅ = 2.971,
    C₆ = 2548.0,
    Pₖ₁ = 1e-12,
    Pₖ₂ = 1e-12,
    Pₖ₃ = 1e-12,
    C₇ = 0.1345,
    C₈ = 351.1,
    C₂₁ = 1e-12,
    C₉ = 0.02869,
    C₁₀ = 1e-12,
    C₂₂ = 1e-12,
    C₁₁ = 0.02928,
    C₁₂ = 4337.0,
    C₂₃ = 1e-12,
    C₁₃ = 0.05098,
    C₁₄ = 476.6,
    C₂₄ = 1e-12,
    C₁₅ = 0.006924,
    C₁₆ = 1e-12,
    C₂₅ = 1e-12,
    C₁₇ = 2.487,
    C₁₈ = 7611.0,
    C₂₆ = 1e-12,
    NK = 2.0,
    ca = 1e-12,
    cb = 1e-12,
    Cx1 = 1.78e6,
    Cx2 = 7.806e3,
    Cdp = 1e-12,
    Cx3 = 5.401e4,
    Cx4 = 8943.0,
    Csp = 1e-12,
    Cx5 = 5.0,
    Cxa = 0.8052,
    Cxb = 3.68,
    Cxc = 4.485,
    Cg1 = 7.41e4,
    Cg2 = 0.8826,
    Cg3 = 1.185e-3,
    a = 1e-12,
    b = 1e-12,
    c = 1e-12,
    pCnuc = 1e-12,
    Tnuc = 1e-12,
    nn = 1e-12,
    Tgrw = 1e-12,
    kr1 = 1e-12,
    krt = 1e-12,
    kr2 = 1e-12,
    kr3 = 1e-12,
    kp1 = 1e-12,
    kpt = 1e-12,
    kp2 = 1e-12,
)

begin
    plt = plot(xlims=(0, 1), ylims=(0, Inf), legendposition=:outerright, widen=1.06)
    for (i, (θ, ψ)) in enumerate(models)
        test = tests[θ]
        # res = ContinuumMechanicsBase.predict(ψ, test, p)
        res = ContinuumMechanicsBase.predict(ψ, test, p; iREXmethod=0, iGSmethod=0)
        # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
        scatter!(plt, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ] ./ 1e6,
                markercolor=i,
                label="$(θ)K:Exp",
            )
        plot!(plt, [first(x) for x in eachcol(res.data.ϵ)], [vonMises(x) for x in eachcol(res.data.σ)],
                linecolor=i,
                label="$(θ)K:Model",
            )
    end
    display(plt)
end

q = ComponentVector(
    C₁ = NaN,
    C₂ = NaN,
    C₃ = p.C₃,
    C₄ = p.C₄,
    C₅ = p.C₅,
    C₆ = p.C₆,
    Pₖ₁ = p.Pₖ₁,
    Pₖ₂ = p.Pₖ₂,
    Pₖ₃ = p.Pₖ₃,
    C₇ = p.C₇,
    C₈ = p.C₈,
    C₂₁ = p.C₂₁,
    C₉ = p.C₉,
    C₁₀ = p.C₁₀,
    C₂₂ = p.C₂₂,
    C₁₁ = p.C₁₁,
    C₁₂ = p.C₁₂,
    C₂₃ = p.C₂₃,
    C₁₃ = p.C₁₃,
    C₁₄ = p.C₁₄,
    C₂₄ = p.C₂₄,
    C₁₅ = p.C₁₅,
    C₁₆ = p.C₁₆,
    C₂₅ = p.C₂₅,
    C₁₇ = p.C₁₇,
    C₁₈ = p.C₁₈,
    C₂₆ = p.C₂₆,
    NK = p.NK,
    ca = p.ca,
    cb = p.cb,
    Cx1 = p.Cx1,
    Cx2 = p.Cx2,
    Cdp = p.Cdp,
    Cx3 = p.Cx3,
    Cx4 = p.Cx4,
    Csp = p.Csp,
    Cx5 = p.Cx5,
    Cxa = p.Cxa,
    Cxb = p.Cxb,
    Cxc = p.Cxc,
    Cg1 = p.Cg1,
    Cg2 = p.Cg2,
    Cg3 = p.Cg3,
    a = p.a,
    b = p.b,
    c = p.c,
    pCnuc = p.pCnuc,
    Tnuc = p.Tnuc,
    nn = p.nn,
    Tgrw = p.Tgrw,
    kr1 = p.kr1,
    krt = p.krt,
    kr2 = p.kr2,
    kr3 = p.kr3,
    kp1 = p.kp1,
    kpt = p.kpt,
    kp2 = p.kp2,
)

begin
    prob = ContinuumMechanicsBase.MaterialOptimizationProblem(
        collect(Cho2019Unified, values(models)),
        collect(BCJMetalUniaxialTest, values(tests)),
        p,
        parameters(first(values(models))),
        AutoFiniteDiff(),
        L2DistLoss();
        ui=q)
    sol = solve(prob, NelderMead()) # 10.044 s (264158926 allocations: 11.01 GiB)
    @btime x=begin # 5.562 ms (165493 allocations: 7.29 MiB) ~ 5.716 ms (165474 allocations: 7.28 MiB)
        pltq = plot(xlims=(0, 1), ylims=(0, Inf), widen=1.06)
        @sync @distributed for (i, (θ, ψ)) in collect(enumerate(models))
            test = tests[θ]
            # calibration = ContinuumMechanicsBase.predict(ψ, test, sol.u)
            calibration = ContinuumMechanicsBase.predict(ψ, test, sol.u; iREXmethod=0, iGSmethod=0)
            # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
            scatter!(pltq, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ] ./ 1e6,
                    markercolor=i,
                    label="$(θ)K:Exp",
                )
            plot!(pltq, [first(x) for x in eachcol(calibration.data.ϵ)], [vonMises(x) for x in eachcol(calibration.data.σ)],
                    linecolor=i,
                    label="$(θ)K:Calib",
                )
        end
    end; display(pltq)

    @show sol.retcode; i, r = 1, deepcopy(q); for (key, value) in zip(keys(p), q)
        if isnan(value)
            r[key] = sol.u[i]
            @printf("\t%s = %.9f,\n", key, r[key])
        end
        global i += 1
    end
end