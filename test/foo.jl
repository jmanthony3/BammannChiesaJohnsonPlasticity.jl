using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, LossFunctions
# using Plots

using Test



df_Tension_e002_295 = CSV.read("Data_Tension_e0002_T295.csv", DataFrame;
        header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
test = BCJMetalUniaxialTest(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, name="exp")
bcj_loading = BCJMetalStrainControl(295.0, 2e-3, float(last(df_Tension_e002_295[!, "Strain"])), 200, :tension)
G = 159e9   # shear modulus [Pa]
μ = 77e9    # bulk modulus [Pa]
ψ = Bammann1990Modeling(bcj_loading, μ)
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
res = ContinuumMechanicsBase.predict(ψ, test, p)
@show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
# plt = scatter(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, label="exp")
# scatter!(plt, [first(x) for x in eachcol(res.data.ϵ)], [symmetricvonMises(x) for x in eachcol(res.data.σ)], label="DK")
# # scatter!(plt, [x[1, 1] for x in res.data.ϵ], [vonMises(x) for x in res.data.σ], label="DK")
# display(plt)

# q = ComponentVector(
#         C₁  = NaN,      C₂  = p.C₂,    # V
#         C₃  = p.C₃,     C₄  = p.C₄,    # Y
#         C₅  = p.C₅,     C₆  = p.C₆,    # f
#         C₇  = p.C₇,     C₈  = p.C₈,    # r_d
#         C₉  = p.C₉,     C₁₀ = p.C₁₀,   # h
#         C₁₁ = p.C₁₁,    C₁₂ = p.C₁₂,   # r_s
#         C₁₃ = p.C₁₃,    C₁₄ = p.C₁₄,   # R_d
#         C₁₅ = p.C₁₅,    C₁₆ = p.C₁₆,   # H
#         C₁₇ = p.C₁₇,    C₁₈ = p.C₁₈,   # R_s
#         C₁₉ = p.C₁₉,    C₂₀ = p.C₂₀    # Y_adj
# )
# prob = BCJPlasticityProblem(ψ, test, p; ad_type=AutoForwardDiff(), ui=q)
# sol = solve(prob, LBFGS())
# # calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
# # scatter!(plt, [first(x) for x in eachcol(calib.data.ϵ)], [symmetricvonMises(x) for x in eachcol(calib.data.σ)], label="DK (Calib.)")
# # # scatter!(plt, [x[1, 1] for x in res.data.ϵ], [vonMises(x) for x in res.data.σ], label="DK")
# # display(plt)
