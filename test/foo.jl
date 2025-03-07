using BammannChiesaJohnsonPlasticity
using ContinuumMechanicsBase
using ComponentArrays
using CSV
using DataFrames
import ForwardDiff
using FiniteDiff
using Optimization, OptimizationOptimJL
using Plots
using Test

params      = begin
    df          = CSV.read("Props_BCJ_4340_fit.csv", DataFrame; header=true, delim=',', types=[String, Float64])
    rowsofconstants = findall(occursin.(r"C\d{2}", df[!, "Comment"]))
    C_0         = Vector{Float64}(undef, 20)
    C_0[rowsofconstants] .= df[!, "For Calibration with vumat"][rowsofconstants]
    bulk_mod    = df[!, "For Calibration with vumat"][findfirst(occursin("Bulk Mod"), df[!, "Comment"])]
    shear_mod   = df[!, "For Calibration with vumat"][findfirst(occursin("Shear Mod"), df[!, "Comment"])]
    ( # collect as dictionary
        C₁  = C_0[ 1],
        C₂  = C_0[ 2],
        C₃  = C_0[ 3],
        C₄  = C_0[ 4],
        C₅  = C_0[ 5],
        C₆  = C_0[ 6],
        C₇  = C_0[ 7],
        C₈  = C_0[ 8],
        C₉  = C_0[ 9],
        C₁₀ = C_0[10],
        C₁₁ = C_0[11],
        C₁₂ = C_0[12],
        C₁₃ = C_0[13],
        C₁₄ = C_0[14],
        C₁₅ = C_0[15],
        C₁₆ = C_0[16],
        C₁₇ = C_0[17],
        C₁₈ = C_0[18],
        C₁₉ = C_0[19],
        C₂₀ = C_0[20],
        G   = bulk_mod,
        μ   = shear_mod)
end
df_Tension_e002_295 = CSV.read("Data_Tension_e0002_T295.csv", DataFrame;
    header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
test = BCJMetalUniaxialTest(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, name="exp")
bcj_loading = BCJMetalStrainControl(295.0, 2e-3, last(df_Tension_e002_295[!, "Strain"]), 200, :tension)
ψ = DK(bcj_loading, params.μ)
p = ComponentVector(
    C₁  = params.C₁,     C₂     = params.C₂,    # V
    C₃  = params.C₃,     C₄     = params.C₄,    # Y
    C₅  = params.C₅,     C₆     = params.C₆,    # f
    C₇  = params.C₇,     C₈     = params.C₈,    # r_d
    C₉  = params.C₉,     C₁₀    = params.C₁₀,   # h
    C₁₁ = params.C₁₁,    C₁₂    = params.C₁₂,   # r_s
    C₁₃ = params.C₁₃,    C₁₄    = params.C₁₄,   # R_d
    C₁₅ = params.C₁₅,    C₁₆    = params.C₁₆,   # H
    C₁₇ = params.C₁₇,    C₁₈    = params.C₁₈,   # R_s
    C₁₉ = params.C₁₉,    C₂₀    = params.C₂₀    # Y_adj
)

# update!(ψ, p)
res = ContinuumMechanicsBase.predict(ψ, test, p)
# [x[1, 1] for x in res.data.λ]
plt = scatter(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, label="exp")
scatter!(plt, [x[1] for x in eachcol(res.data.λ)], [symmetricvonMises(x) for x in eachcol(res.data.s)], label="DK")
# scatter!(plt, [x[1, 1] for x in res.data.λ], [vonMises(x) for x in res.data.s], label="DK")
display(plt)

# q = ComponentVector(
#     C₁  = NaN,          C₂     = 0.99params.C₂,    # V
#     C₃  = NaN,          C₄     = NaN,    # Y
#     C₅  = NaN,          C₆     = NaN,    # f
#     C₇  = NaN,          C₈     = NaN,    # r_d
#     C₉  = NaN,          C₁₀    = NaN,   # h
#     C₁₁ = NaN,          C₁₂    = NaN,   # r_s
#     C₁₃ = NaN,          C₁₄    = NaN,   # R_d
#     C₁₅ = NaN,          C₁₆    = NaN,   # H
#     C₁₇ = NaN,          C₁₈    = NaN    # R_s
# )
# prob = BCJProblem(ψ, test, p; ad_type=AutoFiniteDiff(), ui=q)
# sol = solve(prob, NelderMead())


# # grad=ForwardDiff.gradient(x->sum(ContinuumMechanicsBase.predict(ψ, test, x).data.s[:,50]), p)