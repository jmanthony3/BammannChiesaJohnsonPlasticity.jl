using ContinuumMechanicsBase
using ComponentArrays, StructArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, LossFunctions
using Plots

include("JohnsonCook-functions.jl")

df_Tension_e002_295 = CSV.read("../test/Data_Tension_e0002_T295.csv", DataFrame;
    header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
test = JCUniaxialTest( # (x, y) data from experiment for calibration
    df_Tension_e002_295[!, "Strain"],
    df_Tension_e002_295[!, "Stress"],
    name="exp")
jc_loading = JCStrainControl( # loading conditions
    295.0,                                          # temperature
    2e-3,                                           # strain-rate
    float(last(df_Tension_e002_295[!, "Strain"])),  # final strain
    200)                                            # number of strain increments
Tr = 295.0
Tm = 1793.0
er0 = 1.0
ψ = JC(jc_loading, Tr, Tm, er0)
p = ComponentVector(
    A   = 835.014717457143,
    B   = 810.9199208476191,
    n   = 0.5152083333333339,
    C   = 0.004992526238095374,
    m   = 1.1721130952380956
)

res = ContinuumMechanicsBase.predict(ψ, test, p)
plt = scatter(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"], label="exp")
scatter!(plt, [only(x) for x in eachcol(res.data.ϵ)], [only(x) for x in eachcol(res.data.σ)], label="Johnson-Cook")
q = ComponentVector(
    A   = 835.014717457143,
    B   = 810.9199208476191,
    n   = NaN,
    C   = 0.004992526238095374,
    m   = NaN
)
# prob = MaterialOptimizationProblem(ψ, test, p; ad_type=AutoForwardDiff(), ui=q) # what I had before while testing
prob = ContinuumMechanicsBase.MaterialOptimizationProblem(ψ, test, p, q, AutoForwardDiff(), L2DistLoss()) # with CMB@v0.2.2
sol = solve(prob, LBFGS())
calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
scatter!(plt, [only(x) for x in eachcol(calib.data.ϵ)], [only(x) for x in eachcol(calib.data.σ)], label="JC (Calib.)")