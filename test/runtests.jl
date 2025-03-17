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
    function testmodel(ψ, p, q)
        prob = ContinuumMechanicsBase.MaterialOptimizationProblem(ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(); ui=q)
        return solve(prob, LBFGS())
    end
    @testset "Bammann1990Modeling" begin
        ψ = Bammann1990Modeling(bcj_loading, μ)
        p = ComponentVector(
            C₁  = 3.78679e7,                C₂  = 352.303,      # V
            C₃  = 4.44448e8,                C₄  = 118.644,      # Y
            C₅  = 1.17097,                  C₆  = 644.065,      # f
            C₇  = 1.07958e-8,               C₈  = 5094.6,       # r_s
            C₉  = 9.97698e-10,              C₁₀ = 196.07,       # r_d
            C₁₁ = 7.82595e-10,              C₁₂ = 14102.4,      # R_s
            C₁₃ = 4.37044e-12,              C₁₄ = 7.75843e-12,  # R_d
            C₁₅ = 5.91094e7,                C₁₆ = 0.0022459,    # h
            C₁₇ = 5.27401e8,                C₁₈ = 0.00142715,   # H
        )
        q = ComponentVector(
            C₁  = NaN,      C₂  = NaN,      # V
            C₃  = p.C₃,     C₄  = p.C₄,     # Y
            C₅  = p.C₅,     C₆  = p.C₆,     # f
            C₇  = p.C₇,     C₈  = p.C₈,     # r_s
            C₉  = p.C₉,     C₁₀ = p.C₁₀,    # r_d
            C₁₁ = p.C₁₁,    C₁₂ = p.C₁₂,    # R_s
            C₁₃ = p.C₁₃,    C₁₄ = p.C₁₄,    # R_s
            C₁₅ = p.C₁₅,    C₁₆ = p.C₁₆,    # h
            C₁₇ = p.C₁₇,    C₁₈ = p.C₁₈,    # H
        )
        @test testmodel(ψ, p, q).retcode == SciMLBase.ReturnCode.Success
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
