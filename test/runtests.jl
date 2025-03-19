using BammannChiesaJohnsonPlasticity

using ContinuumMechanicsBase
using ComponentArrays
using CSV, DataFrames
using FiniteDiff
import ForwardDiff
using Optimization, OptimizationOptimJL, LossFunctions

using Test



@testset verbose=true "BammannChiesaJohnsonPlasticity.jl" begin
    df_Tension_e002_295 = CSV.read("Data_Tension_e0002_T295.csv", DataFrame;
        header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
    test = BCJMetalUniaxialTest(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, name="exp")
    bcj_loading = BCJMetalStrainControl(295.0, 2e-3, float(last(df_Tension_e002_295[!, "Strain"])), 200, :tension)
    G = 159e9   # shear modulus [Pa]
    μ = 77e9    # bulk modulus [Pa]
    function testmodel(ψ, test, p, q)
        prob = ContinuumMechanicsBase.MaterialOptimizationProblem(
            ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(); ui=q)
        return solve(prob, LBFGS())
    end
    @testset "Bammann1990Modeling" begin
        ψ = Bammann1990Modeling(bcj_loading, μ)
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
        pred = ContinuumMechanicsBase.predict(ψ, test, p)
        # @show [vonMises(x) for x in eachcol(pred.data.σ)] ./ 1e6
        @test isapprox(31.936, rmse(
            (df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
            ([first(x) for x in eachcol(pred.data.ϵ)], [vonMises(x) for x in eachcol(pred.data.σ)] ./ 1e6)); atol=1e-2)
        q = ComponentVector(
            C₁ = p.C₁,
            C₂ = p.C₂,
            C₃ = p.C₃,
            C₄ = p.C₄,
            C₅ = p.C₅,
            C₆ = p.C₆,
            C₇ = p.C₇,
            C₈ = p.C₈,
            C₉ = p.C₉,
            C₁₀ = p.C₁₀,
            C₁₁ = p.C₁₁,
            C₁₂ = p.C₁₂,
            C₁₃ = p.C₁₃,
            C₁₄ = p.C₁₄,
            C₁₅ = NaN,
            C₁₆ = NaN,
            C₁₇ = NaN,
            C₁₈ = NaN,
        )
        sol = testmodel(ψ, test, p, q)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
        @test isapprox(29.888, rmse(
            (df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
            ([first(x) for x in eachcol(calib.data.ϵ)], [vonMises(x) for x in eachcol(calib.data.σ)] ./ 1e6)); atol=1e-2)
    end

    # # bcj_loading = BCJ_metal(295., 570., 0.15, 200, 1, p)
    # bcj_conf_Tension_e002_295       = kernel(DK, bcj_loading_Tension_e002_295)
    # bcj_ref_Tension_e002_295        = bcj_conf_Tension_e002_295[1]
    # bcj_current_Tension_e002_295    = bcj_conf_Tension_e002_295[2]
    # bcj_history_Tension_e002_295    = bcj_conf_Tension_e002_295[3]
    # σvM = symmetricvonMises(bcj_history_Tension_e002_295.σ__)
    # idx = []
    # for t in df_Tension_e002_295[!, "Strain"]
    #     j = findlast(bcj_history_Tension_e002_295.ϵ__[1, :] .<= t)
    #     push!(idx, !isnothing(j) ? j : findfirst(bcj_history_Tension_e002_295.ϵ__[1, :] .>= t))
    # end
    # err = sum(((df_Tension_e002_295[!, "Stress"] .* 1e6) - σvM[idx]) .^ 2.)
    # x = sum((df_Tension_e002_295[!, "Stress"] .* 1e6) .^ 2)
    # # println(x)
    # err /= if x > 1e9^2
    #     1e9^2
    # elseif x > 1e6^2
    #     1e6^2
    # elseif x > 1e3^2
    #     1e3^2
    # else
    #     1.
    # end
    # # println(err)
    # # p = scatter(df_file[!, "Strain"], df_file[!, "Stress"] .* 1e6, label="Data", ylims=(0., 2e9))
    # # plot!(p, bcj_history.ϵ__[1, :], σvM, label="Model")
    # # display(p)
    # @test isapprox(err, 0.020599315626415465; atol=1e-6)

    # df_Tension_e570_295 = CSV.read("Data_Tension_e570_T295.csv", DataFrame;
    #     header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
    # # bcj_loading = BCJMetalStrainControl(295., 2e-3, 1., 200, 1, p)
    # bcj_loading_Tension_e570_295 = BCJMetalStrainControl(295., 570., 0.15, 200, :tension, p)
    # bcj_conf_Tension_e570_295 = referenceconfiguration(DK, bcj_loading_Tension_e570_295)
    # bcj_ref_Tension_e570_295        = bcj_conf_Tension_e570_295[1]
    # bcj_current_Tension_e570_295    = bcj_conf_Tension_e570_295[2]
    # bcj_history_Tension_e570_295    = bcj_conf_Tension_e570_295[3]
    # solve!(bcj_current_Tension_e570_295, bcj_history_Tension_e570_295)
    # σvM = symmetricvonMises(bcj_history_Tension_e570_295.σ__)
    # idx = []
    # for t in df_Tension_e570_295[!, "Strain"]
    #     j = findlast(bcj_history_Tension_e570_295.ϵ__[1, :] .<= t)
    #     push!(idx, !isnothing(j) ? j : findfirst(bcj_history_Tension_e570_295.ϵ__[1, :] .>= t))
    # end
    # err = sum(((df_Tension_e570_295[!, "Stress"] .* 1e6) - σvM[idx]) .^ 2.)
    # x = sum((df_Tension_e570_295[!, "Stress"] .* 1e6) .^ 2)
    # # println(x)
    # err /= if x > 1e9^2
    #     1e9^2
    # elseif x > 1e6^2
    #     1e6^2
    # elseif x > 1e3^2
    #     1e3^2
    # else
    #     1.
    # end
    # @test isapprox(err, 0.0003468752570447703; atol=1e-6)

    # bcj_conf_comb           = nextloadingphase(
    #     referenceconfiguration(DK, bcj_loading_Tension_e570_295),
    #     bcj_history_Tension_e002_295)
    # # bcj_ref_comb            = bcj_conf_comb[1]
    # # copyto!(bcj_ref_comb, bcj_history_Tension_e002_295)
    # # bcj_current_comb        = bcj_ref_comb
    # # bcj_history_comb        = bcj_conf_comb[3]
    # # record!(bcj_history_comb, 1, bcj_current_comb)
    # bcj_ref_comb        = bcj_conf_comb[1]
    # bcj_current_comb    = bcj_conf_comb[2]
    # bcj_history_comb    = bcj_conf_comb[3]
    # solve!(bcj_current_comb, bcj_history_comb)
    # σvM = symmetricvonMises(bcj_history_comb.σ__)
    # idx = []
    # for t in df_Tension_e570_295[!, "Strain"]
    #     j = findlast(bcj_history_comb.ϵ__[1, :] .<= t)
    #     push!(idx, !isnothing(j) ? j : findfirst(bcj_history_Tension_e570_295.ϵ__[1, :] .>= t))
    # end
    # err = sum(((df_Tension_e570_295[!, "Stress"] .* 1e6) - σvM[idx]) .^ 2.)
    # x = sum((df_Tension_e570_295[!, "Stress"] .* 1e6) .^ 2)
    # # println(x)
    # err /= if x > 1e9^2
    #     1e9^2
    # elseif x > 1e6^2
    #     1e6^2
    # elseif x > 1e3^2
    #     1e3^2
    # else
    #     1.
    # end
    # @test isapprox(err, 2.326137838452948; atol=1e-6)
    # # p = scatter(df_Tension_e570_295[!, "Strain"], df_Tension_e570_295[!, "Stress"] .* 1e6, label="Data", ylims=(0., 2e9))
    # # plot!(p, bcj_history_comb.ϵ__[1, :], σvM, label="Model")
    # # display(p)
    # df_strain = df_Tension_e002_295[!, "Strain"]
    # append!(df_strain, last(df_strain) .+ df_Tension_e570_295[!, "Strain"])
    # df_stress = df_Tension_e002_295[!, "Stress"]
    # append!(df_stress, df_Tension_e570_295[!, "Stress"] .+ (last(df_Tension_e002_295[!, "Stress"]) - first(df_Tension_e570_295[!, "Stress"])))
    # p = scatter(df_strain, df_stress .* 1e6, label="Data", ylims=(0., 2e9))
    # bcj_history_combined = bcj_history_Tension_e002_295 + bcj_history_comb
    # σvM = symmetricvonMises(bcj_history_combined.σ__)
    # idx = []
    # for t in df_strain
    #     j = findlast(bcj_history_combined.ϵ__[1, :] .<= t)
    #     push!(idx, !isnothing(j) ? j : findfirst(bcj_history_combined.ϵ__[1, :] .>= t))
    # end
    # err = sum((df_stress - σvM[idx]) .^ 2.)
    # x = sum(df_stress .^ 2)
    # # println(x)
    # err /= if x > 1e9^2
    #     1e9^2
    # elseif x > 1e6^2
    #     1e6^2
    # elseif x > 1e3^2
    #     1e3^2
    # else
    #     1.
    # end
    # # println(err)
    # plot!(p, bcj_history_combined.ϵ__[1, :], σvM, label="Model")
    # display(p)
    # @test isapprox(err, 3.666486610497823e13; atol=1e6)
end
