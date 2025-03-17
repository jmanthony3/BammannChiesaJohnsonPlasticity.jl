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
ψ = DK(bcj_loading, μ)
p = ComponentVector(
    C₁  = 35016459.896579415,       C₂  = 323.93342698083165,   # V
    C₃  = 500340419.8337271,        C₄  = 143.08381901004486,   # Y
    C₅  = 4.101775377562497,        C₆  = 271.0245526,          # f
    C₇  = 1.0834796217232945e-06,   C₈  = 1023.6278003945317,   # r_s
    C₉  = 2358205093.844017,        C₁₀ = 676421.9935474312,    # h
    C₁₁ = 1.3465080192134937e-10,   C₁₂ = 98.35671405000001,    # r_d
    C₁₃ = 2.533629073577668e-09,    C₁₄ = 403.2291451343492,    # R_s
    C₁₅ = 1159915808.5023918,       C₁₆ = 959557.0948847248,    # H
    C₁₇ = 6.204370386543724e-12,    C₁₈ = 203.95288011132806,   # R_s
    C₁₉ = 1e-10,                    C₂₀ = 1e-10                 # Y_adj
)
res = ContinuumMechanicsBase.predict(ψ, test, p)
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
