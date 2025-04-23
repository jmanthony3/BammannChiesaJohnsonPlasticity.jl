using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, DataInterpolations, LossFunctions

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
prediction = ContinuumMechanicsBase.predict(ψ, test, p)
ϵ = [first(x) for x in test.data.ϵ]
σ = [first(x) for x in test.data.σ]
ϵ̂ = [first(x) for x in eachcol(prediction.data.ϵ)]
σ̂ = [vonMises(x) for x in eachcol(prediction.data.σ)]
# s = linear_interpolation(ϵ̂, σ̂, extrapolation_bc=Line()).(ϵ)
s = CubicSpline(σ̂, ϵ̂; extrapolation=ExtrapolationType.Linear).(ϵ)
# @test isapprox(31.936, rmse(
#     (df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
#     (ϵ, s ./ 1e6)); atol=1e-2)

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
prob = ContinuumMechanicsBase.MaterialOptimizationProblem(
    ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(); ui=q)
sol = solve(prob, LBFGS())

calibration = ContinuumMechanicsBase.predict(ψ, test, sol.u)
ϵ̂ = [first(x) for x in eachcol(calibration.data.ϵ)]
σ̂ = [vonMises(x) for x in eachcol(calibration.data.σ)]
# s = linear_interpolation(ϵ̂, σ̂, extrapolation_bc=Line()).(ϵ)
s = CubicSpline(σ̂, ϵ̂; extrapolation=ExtrapolationType.Linear).(ϵ)
# @test isapprox(29.888, rmse(
#     (df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
#     (ϵ, s ./ 1e6)); atol=1e-2)