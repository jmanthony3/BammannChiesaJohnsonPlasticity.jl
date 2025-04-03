using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, LossFunctions
using Plots

using Test

# include("Cho2019UnifiedStaticDynamic-functions.jl")
include("ext.jl")



df_Fig4a = CSV.read("Cho2019UnifiedStaticDynamic-Fig4a.csv", DataFrame;
        header=true, skipto=3, delim=',', types=Float64)
# test_298K = BCJMetalUniaxialTest( # (x, y) data from experiment for calibration
#     filter(!ismissing, df_Fig4a[!, 1]),
#     filter(!ismissing, df_Fig4a[!, 2]) .* 1e6,
#     name="exp")
# Ω_298K = BCJMetalStrainControl( # loading conditions
#     298.0,                                                      # temperature
#     4e-4,                                                       # strain-rate
#     float(last(filter(!ismissing, df_Fig4a[!, 1]))), # final strain
#     200,                                                        # number of strain increments
#     :tension)                                                   # load direction

ϵ̇ = 4e-4
K = 159e9   # bulk modulus [Pa]
μ = 77e9    # shear modulus [Pa]

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
    models[θ_str] = Bammann1990Modeling(domains[θ_str], μ)
end

tests = sort(tests; rev=false)
domains = sort(domains; rev=false)
models = sort(models; rev=false)

p = ComponentVector(
    C₁ = 9.1402e10,
    C₂ = 258.417,
    C₃ = 1.62805e8,
    C₄ = 363.053,
    C₅ = 1.28544,
    C₆ = 236.047,
    C₇ = 1.04959e-6,
    C₈ = 0.0920373,
    C₉ = 4.07014e-10,
    C₁₀ = 1000.0,
    C₁₁ = 7.07701e-12,
    C₁₂ = 18.6325,
    C₁₃ = 5.07815e-12,
    C₁₄ = 38.7783,
    C₁₅ = 3.77314e7,
    C₁₆ = 0.0111427,
    C₁₇ = 7.87311e6,
    C₁₈ = 0.0155747,
)

# res = ContinuumMechanicsBase.predict(ψ, test, p)
# @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
RES = []
begin
    plt = plot(xlims=(0, 1), ylims=(0, Inf), widen=1.06)
    for (i, (θ, ψ)) in enumerate(models)
        test = tests[θ]
        res = ContinuumMechanicsBase.predict(ψ, test, p)
        push!(RES, res)
        # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
        scatter!(plt, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ],
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
        C₁  = NaN,      C₂  = NaN,    # V
        C₃  = NaN,      C₄  = NaN,    # Y
        C₅  = NaN,      C₆  = NaN,    # f
        C₇  = p.C₇,     C₈  = p.C₈,    # r_d
        C₉  = p.C₉,     C₁₀ = p.C₁₀,   # h
        C₁₁ = p.C₁₁,    C₁₂ = p.C₁₂,   # r_s
        C₁₃ = p.C₁₃,    C₁₄ = p.C₁₄,   # R_d
        C₁₅ = p.C₁₅,    C₁₆ = p.C₁₆,   # H
        C₁₇ = p.C₁₇,    C₁₈ = p.C₁₈,   # R_s
)
# sol = Dict()
# for (i, (θ, ψ)) in enumerate(models)
#     test = tests[θ]
#     prob = ContinuumMechanicsBase.MaterialOptimizationProblem(ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(), ui=q)
#     sol[θ] = solve(prob, LBFGS())
# end
prob = ContinuumMechanicsBase.MaterialOptimizationProblem(
    collect(Bammann1990Modeling, values(models)),
    collect(BCJMetalUniaxialTest, values(tests)),
    p,
    parameters(first(values(models))),
    AutoForwardDiff(),
    L2DistLoss();
    ui=q)
# prob = BCJPlasticityProblem(ψ, test, p; ad_type=AutoForwardDiff(), ui=q)
sol = solve(prob, LBFGS())
# # calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
# # scatter!(plt, [first(x) for x in eachcol(calib.data.ϵ)], [symmetricvonMises(x) for x in eachcol(calib.data.σ)], label="DK (Calib.)")
# # # scatter!(plt, [x[1, 1] for x in res.data.ϵ], [vonMises(x) for x in res.data.σ], label="DK")
# # display(plt)
