using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, Interpolations, LossFunctions

using Test

function MMO(
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
        ŝ = σ # vonMises.(σ̂) # linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵ̂)
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



df_Tension_e002_295 = CSV.read("Data_Tension_e0002_T295.csv", DataFrame;
    header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
test = BCJMetalUniaxialTest(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, name="exp")
Ω = BCJMetalStrainControl(295.0, 2e-3, float(last(df_Tension_e002_295[!, "Strain"])), 200, :tension)
K = 159e9   # bulk modulus [Pa]
μ = 77e9    # shear modulus [Pa]
ψ = Bammann1990Modeling(Ω, μ)
p = ComponentVector(
    C₁ = 9.98748e10,
    C₂ = 1483.14,
    C₃ = 1.61687e8,
    C₄ = 382.443,
    C₅ = 1.65237,
    C₆ = 1320.97,
    C₇ = 0.000195306,
    C₈ = 1504.62,
    C₉ = 4.04209e-10,
    C₁₀ = 993.109,
    C₁₁ = 7.02824e-12,
    C₁₂ = 18.5041,
    C₁₃ = 5.04316e-9,
    C₁₄ = 2153.13,
    C₁₅ = 3.73042e7,
    C₁₆ = 1792.72,
    C₁₇ = 9.56827e6,
    C₁₈ = 1214.34,
)

q = ComponentVector(
    C₁ = p.C₁,      C₂ = p.C₂,
    C₃ = p.C₃,      C₄ = p.C₄,
    C₅ = p.C₅,      C₆ = p.C₆,
    C₇ = p.C₇,      C₈ = p.C₈,
    C₉ = p.C₉,      C₁₀ = p.C₁₀,
    C₁₁ = p.C₁₁,    C₁₂ = p.C₁₂,
    C₁₃ = p.C₁₃,    C₁₄ = p.C₁₄,
    C₁₅ = NaN,      C₁₆ = NaN,
    C₁₇ = NaN,      C₁₈ = NaN,
)

prob = MMO(ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(); ui=q)
sol = solve(prob, LBFGS())

# sol = testmodel(ψ, test, p, q)
# @test sol.retcode == SciMLBase.ReturnCode.Success
# calibration = ContinuumMechanicsBase.predict(ψ, test, sol.u)
# ϵ̂ = [first(x) for x in eachcol(calibration.data.ϵ)]
# σ̂ = [vonMises(x) for x in eachcol(calibration.data.σ)]
# s = linear_interpolation(ϵ̂, σ̂, extrapolation_bc=Line()).(ϵ)
# @test isapprox(29.888, rmse(
#     (df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
#     (ϵ, s ./ 1e6)); atol=1e-2)