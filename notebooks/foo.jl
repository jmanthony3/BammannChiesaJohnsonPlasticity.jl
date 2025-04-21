using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, LossFunctions
using Plots
using Interpolations

using Test

include("Cho2019UnifiedStaticDynamic-functions.jl")
include("ext.jl")



df_Tension_e002_295 = CSV.read("../test/Data_Tension_e0002_T295.csv", DataFrame;
    header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
test = BCJMetalUniaxialTest(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, name="exp")
Ω = BCJMetalStrainControl(295.0, 2e-3, float(last(df_Tension_e002_295[!, "Strain"])), 200, :tension)
K = 159e9   # bulk modulus [Pa]
μ = 77e9    # shear modulus [Pa]

ψ = Bammann1990Modeling(Ω, μ)

p0 = ComponentVector(
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
prediction = ContinuumMechanicsBase.predict(ψ, test, p0)

ϵ = [first(x) for x in test.data.ϵ]
σ = [first(x) for x in test.data.σ]
ϵ̂ = [first(x) for x in eachcol(prediction.data.ϵ)]
σ̂ = collect(eachcol(prediction.data.σ))
# resϵ = [x[1, 1] for x in pred.data.ϵ]
# testϵ = [x[1, 1] for x in test.data.ϵ]
# s = collect([[x...] for x in eachcol(pred.data.σ)[[findlast(x .>= resϵ) for x in testϵ]]])
# # s = collect([[x...] for x in pred.data.σ[[findlast(x .>= resϵ) for x in testϵ]]])
ŝ = linear_interpolation(ϵ, σ, extrapolation_bc=Line()).(ϵ̂)
err = map(i -> loss.(i[1], vonMises(i[2])), zip(ŝ, σ̂)) |> mean

# # res = ContinuumMechanicsBase.predict(ψ, test, p)
# # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
# RES = []
# begin
#     plt = plot(xlims=(0, 1), ylims=(0, Inf), widen=1.06)
#     for (i, (θ, ψ)) in enumerate(models)
#         test = tests[θ]
#         res = ContinuumMechanicsBase.predict(ψ, test, p)
#         push!(RES, res)
#         # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
#         scatter!(plt, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ],
#                 markercolor=i,
#                 label="$(θ)K:Exp",
#             )
#         plot!(plt, [first(x) for x in eachcol(res.data.ϵ)], [vonMises(x) for x in eachcol(res.data.σ)],
#                 linecolor=i,
#                 label="$(θ)K:Model",
#             )
#     end
#     display(plt)
# end

# q = ComponentVector(
#         C₁  = NaN,      C₂  = NaN,    # V
#         C₃  = NaN,      C₄  = NaN,    # Y
#         C₅  = NaN,      C₆  = NaN,    # f
#         C₇  = p.C₇,     C₈  = p.C₈,    # r_d
#         C₉  = p.C₉,     C₁₀ = p.C₁₀,   # h
#         C₁₁ = p.C₁₁,    C₁₂ = p.C₁₂,   # r_s
#         C₁₃ = p.C₁₃,    C₁₄ = p.C₁₄,   # R_d
#         C₁₅ = p.C₁₅,    C₁₆ = p.C₁₆,   # H
#         C₁₇ = p.C₁₇,    C₁₈ = p.C₁₈,   # R_s
# )
# # sol = Dict()
# # for (i, (θ, ψ)) in enumerate(models)
# #     test = tests[θ]
# #     prob = ContinuumMechanicsBase.MaterialOptimizationProblem(ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(), ui=q)
# #     sol[θ] = solve(prob, LBFGS())
# # end
# prob = ContinuumMechanicsBase.MaterialOptimizationProblem(
#     collect(Bammann1990Modeling, values(models)),
#     collect(BCJMetalUniaxialTest, values(tests)),
#     p,
#     parameters(first(values(models))),
#     AutoForwardDiff(),
#     L2DistLoss();
#     ui=q)
# # prob = BCJPlasticityProblem(ψ, test, p; ad_type=AutoForwardDiff(), ui=q)
# sol = solve(prob, LBFGS())
# # # calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
# # # scatter!(plt, [first(x) for x in eachcol(calib.data.ϵ)], [symmetricvonMises(x) for x in eachcol(calib.data.σ)], label="DK (Calib.)")
# # # # scatter!(plt, [x[1, 1] for x in res.data.ϵ], [vonMises(x) for x in res.data.σ], label="DK")
# # # display(plt)
