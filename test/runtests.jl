using BammannChiesaJohnsonPlasticity
using CSV
using DataFrames
using Plots
using Test

@testset verbose=true "BammannChiesaJohnsonPlasticity.jl" begin
    params      = begin
        df          = CSV.read("Props_BCJ_4340_fit.csv", DataFrame; header=true, delim=',', types=[String, Float64])
        rowsofconstants = findall(occursin.(r"C\d{2}", df[!, "Comment"]))
        C_0         = Vector{Float64}(undef, 20)
        C_0[rowsofconstants] .= df[!, "For Calibration with vumat"][rowsofconstants]
        bulk_mod    = df[!, "For Calibration with vumat"][findfirst(occursin("Bulk Mod"), df[!, "Comment"])]
        shear_mod   = df[!, "For Calibration with vumat"][findfirst(occursin("Shear Mod"), df[!, "Comment"])]
        Dict( # collect as dictionary
            "C01"       => C_0[ 1],
            "C02"       => C_0[ 2],
            "C03"       => C_0[ 3],
            "C04"       => C_0[ 4],
            "C05"       => C_0[ 5],
            "C06"       => C_0[ 6],
            "C07"       => C_0[ 7],
            "C08"       => C_0[ 8],
            "C09"       => C_0[ 9],
            "C10"       => C_0[10],
            "C11"       => C_0[11],
            "C12"       => C_0[12],
            "C13"       => C_0[13],
            "C14"       => C_0[14],
            "C15"       => C_0[15],
            "C16"       => C_0[16],
            "C17"       => C_0[17],
            "C18"       => C_0[18],
            "C19"       => C_0[19],
            "C20"       => C_0[20],
            "bulk_mod"  => bulk_mod,
            "shear_mod" => shear_mod)
    end
    df_Tension_e002_295 = CSV.read("Data_Tension_e0002_T295.csv", DataFrame;
        header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
    bcj_loading_Tension_e002_295 = BCJMetalStrainControl(295., 2e-3, 1., 200, 1, params)
    # bcj_loading = BCJ_metal(295., 570., 0.15, 200, 1, params)
    bcj_conf_Tension_e002_295 = referenceconfiguration(DK, bcj_loading_Tension_e002_295)
    bcj_ref_Tension_e002_295        = bcj_conf_Tension_e002_295[1]
    bcj_current_Tension_e002_295    = bcj_conf_Tension_e002_295[2]
    bcj_history_Tension_e002_295    = bcj_conf_Tension_e002_295[3]
    solve!(bcj_current_Tension_e002_295, bcj_history_Tension_e002_295)
    σvM = symmetricvonMises(bcj_history_Tension_e002_295.σ__)
    idx = []
    for t in df_Tension_e002_295[!, "Strain"]
        j = findlast(bcj_history_Tension_e002_295.ϵ__[1, :] .<= t)
        push!(idx, !isnothing(j) ? j : findfirst(bcj_history_Tension_e002_295.ϵ__[1, :] .>= t))
    end
    err = sum(((df_Tension_e002_295[!, "Stress"] .* 1e6) - σvM[idx]) .^ 2.)
    x = sum((df_Tension_e002_295[!, "Stress"] .* 1e6) .^ 2)
    # println(x)
    err /= if x > 1e9^2
        1e9^2
    elseif x > 1e6^2
        1e6^2
    elseif x > 1e3^2
        1e3^2
    else
        1.
    end
    # println(err)
    # p = scatter(df_file[!, "Strain"], df_file[!, "Stress"] .* 1e6, label="Data", ylims=(0., 2e9))
    # plot!(p, bcj_history.ϵ__[1, :], σvM, label="Model")
    # display(p)
    @test isapprox(err, 0.020599315626415465; atol=1e-6)

    df_Tension_e570_295 = CSV.read("Data_Tension_e570_T295.csv", DataFrame;
        header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
    # bcj_loading = BCJMetalStrainControl(295., 2e-3, 1., 200, 1, params)
    bcj_loading_Tension_e570_295 = BCJMetalStrainControl(295., 570., 0.15, 200, 1, params)
    bcj_conf_Tension_e570_295 = referenceconfiguration(DK, bcj_loading_Tension_e570_295)
    bcj_ref_Tension_e570_295        = bcj_conf_Tension_e570_295[1]
    bcj_current_Tension_e570_295    = bcj_conf_Tension_e570_295[2]
    bcj_history_Tension_e570_295    = bcj_conf_Tension_e570_295[3]
    solve!(bcj_current_Tension_e570_295, bcj_history_Tension_e570_295)
    σvM = symmetricvonMises(bcj_history_Tension_e570_295.σ__)
    idx = []
    for t in df_Tension_e570_295[!, "Strain"]
        j = findlast(bcj_history_Tension_e570_295.ϵ__[1, :] .<= t)
        push!(idx, !isnothing(j) ? j : findfirst(bcj_history_Tension_e570_295.ϵ__[1, :] .>= t))
    end
    err = sum(((df_Tension_e570_295[!, "Stress"] .* 1e6) - σvM[idx]) .^ 2.)
    x = sum((df_Tension_e570_295[!, "Stress"] .* 1e6) .^ 2)
    # println(x)
    err /= if x > 1e9^2
        1e9^2
    elseif x > 1e6^2
        1e6^2
    elseif x > 1e3^2
        1e3^2
    else
        1.
    end
    @test isapprox(err, 0.0003468752570447703; atol=1e-6)

    bcj_conf_comb           = referenceconfiguration(DK, bcj_loading_Tension_e570_295)
    bcj_ref_comb            = bcj_conf_comb[1]
    copyto!(bcj_ref_comb, bcj_history_Tension_e002_295)
    bcj_current_comb        = bcj_ref_comb
    bcj_history_comb        = bcj_conf_comb[3]
    record!(bcj_history_comb, 1, bcj_current_comb)
    solve!(bcj_current_comb, bcj_history_comb)
    σvM = symmetricvonMises(bcj_history_comb.σ__)
    idx = []
    for t in df_Tension_e570_295[!, "Strain"]
        j = findlast(bcj_history_comb.ϵ__[1, :] .<= t)
        push!(idx, !isnothing(j) ? j : findfirst(bcj_history_Tension_e570_295.ϵ__[1, :] .>= t))
    end
    err = sum(((df_Tension_e570_295[!, "Stress"] .* 1e6) - σvM[idx]) .^ 2.)
    x = sum((df_Tension_e570_295[!, "Stress"] .* 1e6) .^ 2)
    # println(x)
    err /= if x > 1e9^2
        1e9^2
    elseif x > 1e6^2
        1e6^2
    elseif x > 1e3^2
        1e3^2
    else
        1.
    end
    @test isapprox(err, 2.326137838452948; atol=1e-6)
    # p = scatter(df_Tension_e570_295[!, "Strain"], df_Tension_e570_295[!, "Stress"] .* 1e6, label="Data", ylims=(0., 2e9))
    # plot!(p, bcj_history_comb.ϵ__[1, :], σvM, label="Model")
    # display(p)
    df_strain = df_Tension_e002_295[!, "Strain"]
    append!(df_strain, last(df_strain) .+ df_Tension_e570_295[!, "Strain"])
    df_stress = df_Tension_e002_295[!, "Stress"]
    append!(df_stress, df_Tension_e570_295[!, "Stress"] .+ (last(df_Tension_e002_295[!, "Stress"]) - first(df_Tension_e570_295[!, "Stress"])))
    p = scatter(df_strain, df_stress .* 1e6, label="Data", ylims=(0., 2e9))
    bcj_history_combined = bcj_history_Tension_e002_295 + bcj_history_comb
    σvM = symmetricvonMises(bcj_history_combined.σ__)
    idx = []
    for t in df_strain
        j = findlast(bcj_history_combined.ϵ__[1, :] .<= t)
        push!(idx, !isnothing(j) ? j : findfirst(bcj_history_combined.ϵ__[1, :] .>= t))
    end
    err = sum((df_stress - σvM[idx]) .^ 2.)
    x = sum(df_stress .^ 2)
    # println(x)
    err /= if x > 1e9^2
        1e9^2
    elseif x > 1e6^2
        1e6^2
    elseif x > 1e3^2
        1e3^2
    else
        1.
    end
    # println(err)
    plot!(p, bcj_history_combined.ϵ__[1, :], σvM, label="Model")
    display(p)
    @test isapprox(err, 3.666486610497823e13; atol=1e6)
end
