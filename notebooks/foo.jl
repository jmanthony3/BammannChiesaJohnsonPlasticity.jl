using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, Interpolations, LossFunctions
using Plots

include("Cho2019UnifiedStaticDynamic-functions.jl")

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
    Pₖ₁ = 0.0,
    Pₖ₂ = 0.0,
    Pₖ₃ = 0.0,
    C₇ = 0.1345,
    C₈ = 351.1,
    C₂₁ = 0.0,
    C₉ = 0.02869,
    C₁₀ = 0.0,
    C₂₂ = 0.0,
    C₁₁ = 0.02928,
    C₁₂ = 4337.0,
    C₂₃ = 0.0,
    C₁₃ = 0.05098,
    C₁₄ = 476.6,
    C₂₄ = 0.0,
    C₁₅ = 0.006924,
    C₁₆ = 0.0,
    C₂₅ = 0.0,
    C₁₇ = 2.487,
    C₁₈ = 7611.0,
    C₂₆ = 0.0,
    NK = 2.0,
    ca = 0.0,
    cb = 0.0,
    Cx1 = 1.78e6,
    Cx2 = 7.806e3,
    Cdp = 0.0,
    Cx3 = 5.401e4,
    Cx4 = 8943.0,
    Csp = 0.0,
    Cx5 = 5.0,
    Cxa = 0.8052,
    Cxb = 3.68,
    Cxc = 4.485,
    Cg1 = 7.41e4,
    Cg2 = 0.8826,
    Cg3 = 1.185e-3,
    a = 0.0,
    b = 0.0,
    c = 0.0,
    pCnuc = 0.0,
    Tnuc = 0.0,
    nn = 0.0,
    Tgrw = 0.0,
    kr1 = 0.0,
    krt = 0.0,
    kr2 = 0.0,
    kr3 = 0.0,
    kp1 = 0.0,
    kpt = 0.0,
    kp2 = 0.0,
)

begin
    plt = plot(xlims=(0, 1), ylims=(0, Inf), legendposition=:outerright, widen=1.06)
    for (i, (θ, ψ)) in enumerate(models)
        test = tests[θ]
        res = ContinuumMechanicsBase.predict(ψ, test, p)
        # res = ContinuumMechanicsBase.predict(ψ, test, p; iREXmethod=0, iGSmethod=0)
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

# pltq = plot(xlims=(0, 1), ylims=(0, Inf), widen=1.06)

# q = ComponentVector(
#     C₁ = NaN,
#     C₂ = NaN,
#     C₃ = p.C₃,
#     C₄ = p.C₄,
#     C₅ = p.C₅,
#     C₆ = p.C₆,
#     Pₖ₁ = p.Pₖ₁,
#     Pₖ₂ = p.Pₖ₂,
#     Pₖ₃ = p.Pₖ₃,
#     C₇ = p.C₇,
#     C₈ = p.C₈,
#     C₂₁ = p.C₂₁,
#     C₉ = p.C₉,
#     C₁₀ = p.C₁₀,
#     C₂₂ = p.C₂₂,
#     C₁₁ = p.C₁₁,
#     C₁₂ = p.C₁₂,
#     C₂₃ = p.C₂₃,
#     C₁₃ = p.C₁₃,
#     C₁₄ = p.C₁₄,
#     C₂₄ = p.C₂₄,
#     C₁₅ = p.C₁₅,
#     C₁₆ = p.C₁₆,
#     C₂₅ = p.C₂₅,
#     C₁₇ = p.C₁₇,
#     C₁₈ = p.C₁₈,
#     C₂₆ = p.C₂₆,
#     NK = p.NK,
#     ca = p.ca,
#     cb = p.cb,
#     Cx1 = p.Cx1,
#     Cx2 = p.Cx2,
#     Cdp = p.Cdp,
#     Cx3 = p.Cx3,
#     Cx4 = p.Cx4,
#     Csp = p.Csp,
#     Cx5 = p.Cx5,
#     Cxa = p.Cxa,
#     Cxb = p.Cxb,
#     Cxc = p.Cxc,
#     n = p.n,
#     ω₀ = p.ω₀,
#     Cg1 = p.Cg1,
#     Cg2 = p.Cg2,
#     Cg3 = p.Cg3,
#     z = p.z,
#     a = p.a,
#     b = p.b,
#     c = p.c,
#     pCnuc = p.pCnuc,
#     Tnuc = p.Tnuc,
#     nn = p.nn,
#     Tgrw = p.Tgrw,
#     kr1 = p.kr1,
#     krt = p.krt,
#     kr2 = p.kr2,
#     kr3 = p.kr3,
#     kp1 = p.kp1,
#     kpt = p.kpt,
#     kp2 = p.kp2,
# )

# prob = ContinuumMechanicsBase.MaterialOptimizationProblem(
#     collect(Cho2019Unified, values(models)),
#     collect(BCJMetalUniaxialTest, values(tests)),
#     p,
#     parameters(first(values(models))),
#     AutoForwardDiff(),
#     L2DistLoss();
#     ui=q)
# sol = solve(prob, LBFGS())
# for (i, (θ, ψ)) in enumerate(models)
#     test = tests[θ]
#     calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
#     # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
#     scatter!(pltq, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ],
#             markercolor=i,
#             label="$(θ)K:Exp",
#         )
#     plot!(pltq, [first(x) for x in eachcol(calib.data.ϵ)], [vonMises(x) for x in eachcol(calib.data.σ)],
#             linecolor=i,
#             label="$(θ)K:Calib",
#         )
# end

# @show sol.retcode; i, r = 1, deepcopy(q); for (key, value) in zip(keys(p), q)
#     if isnan(value)
#         r[key] = sol.u[i]
#         @printf("\t%s = %.9f,\n", key, r[key])
#     end
#     global i += 1
# end; @show r