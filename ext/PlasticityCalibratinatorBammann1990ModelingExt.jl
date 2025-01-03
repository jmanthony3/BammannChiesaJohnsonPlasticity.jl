using CSV
using DataFrames
using Distributed
using GLMakie
using LaTeXStrings

set_theme!(theme_latexfonts())

materialprops_Bammann1990Modeling = Dict(
    # Comment,For Calibration with vumat
    "C01"       => 4.959788864230217e7,
    "C02"       => 442.6117395684036,
    "C03"       => 4.532613412908138e8,
    "C04"       => 110.28482546319775,
    "C05"       => 7.452174484760692,
    "C06"       => 647.0405985495495,
    "C07"       => 1.0845641859092038e-8,
    "C08"       => 5118.139001972659,
    "C09"       => 1.546170506974806e9,
    "C10"       => 0.0,
    "C11"       => 0.0,
    "C12"       => 0.0,
    "C13"       => 7.862112240331102e-10,
    "C14"       => 1416.7510504720378,
    "C15"       => 4.876522918628674e8,
    "C16"       => 0.0,
    "C17"       => 7.794279114226601e-12,
    "C18"       => 291.9445631223215,
    "Bulk Mod"  => 159000000000.0,
    "Shear Mod" => 77000000000.0,
)
materialconsts_Bammann1990Modeling = collect(
    "C01"       => 4.959788864230217e7,
    "C02"       => 442.6117395684036,
    "C03"       => 4.532613412908138e8,
    "C04"       => 110.28482546319775,
    "C05"       => 7.452174484760692,
    "C06"       => 647.0405985495495,
    "C07"       => 1.0845641859092038e-8,
    "C08"       => 5118.139001972659,
    "C09"       => 1.546170506974806e9,
    "C10"       => 0.0,
    "C11"       => 0.0,
    "C12"       => 0.0,
    "C13"       => 7.862112240331102e-10,
    "C14"       => 1416.7510504720378,
    "C15"       => 4.876522918628674e8,
    "C16"       => 0.0,
    "C17"       => 7.794279114226601e-12,
    "C18"       => 291.9445631223215,
)

PlasticityCalibratinator.materialproperties(::Type{Bammann1990Modeling}) = materialprops_Bammann1990Modeling
PlasticityCalibratinator.materialconstants(::Type{Bammann1990Modeling}) = materialconsts_Bammann1990Modeling

PlasticityCalibratinator.characteristicequations(::Type{Bammann1990Modeling}) = [
    # plasticstrainrate
    L"\dot{\epsilon}_{p} = f(\theta)\sinh\left[ \frac{ \{|\mathbf{\xi}| - \kappa - Y(\theta) \} }{ V(\theta) } \right]\frac{\mathbf{\xi}'}{|\mathbf{\xi}'|}\text{, let }\mathbf{\xi}' = \mathbf{\sigma}' - \mathbf{\alpha}'"
    # kinematichardening
    L"\dot{\mathbf{\alpha}} = h\mu(\theta)\dot{\epsilon}_{p} - [r_{d}(\theta)|\dot{\epsilon}_{p}| + r_{s}(\theta)]|\mathbf{\alpha}|\mathbf{\alpha}"
    # isotropichardening
    L"\dot{\kappa} = H\mu(\theta)\dot{\epsilon}_{p} - [R_{d}(\theta)|\dot{\epsilon}_{p}| + R_{s}(\theta)]\kappa^{2}"
    # flowrule
    L"F = |\sigma - \alpha| - \kappa - \beta(|\dot{\epsilon}_{p}|, \theta)"
    # initialyieldstressbeta
    L"\beta(\dot{\epsilon}_{p}, \theta) = Y(\theta) + V(\theta)\sinh^{-1}\left(\frac{|\dot{\epsilon}_{p}|}{f(\theta)}\right)"
]
PlasticityCalibratinator.dependenceequations(::Type{Bammann1990Modeling})     = [
    L"V = C_{ 1} \mathrm{exp}(-C_{ 2} / \theta)",       # V
    L"Y = C_{ 3} \mathrm{exp}( C_{ 4} / \theta)",       # Y
    L"f = C_{ 5} \mathrm{exp}(-C_{ 6} / \theta)",       # f
    L"r_{d} = C_{ 7} \mathrm{exp}(-C_{ 8} / \theta)",   # rd
    L"h = C_{ 9} \mathrm{exp}( C_{10} / \theta)",       # h
    L"r_{s} = C_{11} \mathrm{exp}(-C_{12} / \theta)",   # rs
    L"R_{d} = C_{13} \mathrm{exp}(-C_{14} / \theta)",   # Rd
    L"H = C_{15} \mathrm{exp}( C_{16} / \theta)",       # H
    L"R_{s} = C_{17} \mathrm{exp}(-C_{18} / \theta)",   # Rs
]
PlasticityCalibratinator.dependencesliders(::Type{Bammann1990Modeling})       = [
    # V
    [
        (label=L"C_{ 1}", range=range(0., 5materialconsts_Bammann1990Modeling["C01"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C01"]),
        (label=L"C_{ 2}", range=range(0., 5materialconsts_Bammann1990Modeling["C02"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C02"]),
    ],
    # Y
    [
        (label=L"C_{ 3}", range=range(0., 5materialconsts_Bammann1990Modeling["C03"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C03"]),
        (label=L"C_{ 4}", range=range(0., 5materialconsts_Bammann1990Modeling["C04"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C04"]),
    ],
    # f
    [
        (label=L"C_{ 5}", range=range(0., 5materialconsts_Bammann1990Modeling["C05"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C05"]),
        (label=L"C_{ 6}", range=range(0., 5materialconsts_Bammann1990Modeling["C06"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C06"]),
    ],
    # rd
    [
        (label=L"C_{ 7}", range=range(0., 5materialconsts_Bammann1990Modeling["C07"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C07"]),
        (label=L"C_{ 8}", range=range(0., 5materialconsts_Bammann1990Modeling["C08"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C08"]),
    ],
    # h
    [
        (label=L"C_{ 9}", range=range(0., 5materialconsts_Bammann1990Modeling["C09"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C09"]),
        (label=L"C_{10}", range=range(0., 5materialconsts_Bammann1990Modeling["C10"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C10"]),
    ],
    # rs
    [
        (label=L"C_{11}", range=range(0., 5materialconsts_Bammann1990Modeling["C11"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C11"]),
        (label=L"C_{12}", range=range(0., 5materialconsts_Bammann1990Modeling["C12"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C12"]),
    ],
    # Rd
    [
        (label=L"C_{13}", range=range(0., 5materialconsts_Bammann1990Modeling["C13"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C13"]),
        (label=L"C_{14}", range=range(0., 5materialconsts_Bammann1990Modeling["C14"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C14"]),
    ],
    # H
    [
        (label=L"C_{15}", range=range(0., 5materialconsts_Bammann1990Modeling["C15"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C15"]),
        (label=L"C_{16}", range=range(0., 5materialconsts_Bammann1990Modeling["C16"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C16"]),
    ],
    # Rs
    [
        (label=L"C_{17}", range=range(0., 5materialconsts_Bammann1990Modeling["C17"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C17"]),
        (label=L"C_{18}", range=range(0., 5materialconsts_Bammann1990Modeling["C18"]; length=1_000), format="{:.3e}", startvalue=materialconsts_Bammann1990Modeling["C18"]),
    ],
]

################################################################
#  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  #
# PlasticityBase
#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  #
# PlasticityCalibratinator
#  v  v  v  v  v  v  v  v  v  v  v  v  v  v  v  v  v  v  v  v  #
################################################################

# FORMATTING: test_data[i][[e_data,s_data] , [e_model,s_model] , [e_err, s_err]]
function PlasticityCalibratinator.modeldata(::Type{Bammann1990Modeling}, modelinputs::ModelInputs, params)::ModelData
    # files, incnum, params, Scale_MPa
    files = modelinputs.expdatasets
    incnum = modelinputs.incnum
    Scale_MPa = modelinputs.stressscale
    # istate = Int64(modelinputs.loading_axial)
    loadstate = modelinputs.loading_axial ? :tension : :compression
    test_cond   = Dict(
        "StrainRate"    => [],
        "Temp"          => [],
        "Name"          => []
    )
    test_data   = Dict(
        "Data_E"        => [],
        "Data_S"        => [],
        "Model_E"       => [],
        "Model_S"       => [],
        "Model_alph"    => [],
        "Model_kap"     => [],
        "Model_VM"      => []
    )
    nsets = 0
    if first(files) != ""
        nsets = length(files)
        for (i, file) in enumerate(files)
            df_file = CSV.read(file, DataFrame; header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])

            # add stress/strain data:
            strn = df_file[!, "Strain"]
            strs = df_file[!, "Stress"] .* Scale_MPa
            er   = df_file[!, "Strain Rate"]
            T    = df_file[!, "Temperature"]
            name = df_file[!, "Name"]
            # check data entered
            if length(strn) != length(strs)
                println("ERROR! Data from  '", file , "'  has bad stress-strain data lengths")
            end

            # store the stress-strain data
            push!(test_cond["StrainRate"],  (i, first(er)))
            push!(test_cond["Temp"],        (i, first(T)))
            push!(test_cond["Name"],        (i, first(name)))
            push!(test_data["Data_E"],      (i, strn))
            push!(test_data["Data_S"],      (i, strs))
        end
        for key in ("StrainRate", "Temp", "Name")
            for (i, x) in enumerate(sort(test_cond[key], by=x->first(x)))
                test_cond[key][i] = last(x)
            end
        end
        for key in ("Data_E", "Data_S")
            for (i, x) in enumerate(sort(test_data[key], by=x->first(x)))
                test_data[key][i] = last(x)
            end
        end
        @sync @distributed for i in range(1, nsets)
        # for i in range(1, nsets)
            # test_data, test_cond, incnum, istate, params, i, ISV_Model
            sol = plotdata_updatekernel(Bammann1990Modeling, test_data, test_cond, incnum, loadstate, params, i)
            # test_data[i][1] = [E,S,Al,kap,tot,SVM]             #Store model stress/strain data
            # push!(test_data["Model_E"],     E)
            # push!(test_data["Model_S"],     S)
            # push!(test_data["Model_alph"],  Al)
            # push!(test_data["Model_kap"],   κ)
            # push!(test_data["Model_tot"],   Tot)
            # push!(test_data["Model_VM"],    SVM)
            push!(test_data["Model_E"],     (i, sol.ϵ))
            push!(test_data["Model_S"],     (i, sol.σ))
            push!(test_data["Model_alph"],  (i, sol.α))
            push!(test_data["Model_kap"],   (i, sol.κ))
            push!(test_data["Model_VM"],    (i, sol.σvM))
        end
        for key in ("Model_E", "Model_S", "Model_alph", "Model_kap", "Model_VM")
            for (i, x) in enumerate(sort(test_data[key], by=x->first(x)))
                test_data[key][i] = last(x)
            end
        end
    end
    return ModelData{Bammann1990Modeling}(Bammann1990Modeling, modelinputs, nsets, test_data, test_cond, params, deepcopy(materialconstants(Bammann1990Modeling)), deepcopy(materialconstants(Bammann1990Modeling)), incnum, Scale_MPa)
end

function PlasticityCalibratinator.plotdata_initialize(::Type{Bammann1990Modeling}, nsets, test_data)
    plot_isvs = []
    dataseries = if !isempty(plot_isvs)
        [[], [], [], []]
    else
        [[], []]
    end
    if nsets > 0
        for i in range(1, nsets)
            push!(dataseries[1], Observable(DataFrame(x=test_data["Data_E"][i], y=test_data["Data_S"][i])))
            push!(dataseries[2], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_VM"][i])))
            if !isempty(plot_isvs)
                push!(dataseries[3], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_alph"][i])))
                push!(dataseries[4], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_kap"][i])))
                # push!(dataseries[5], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_tot"][i])))
                # push!(dataseries[6], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_S"][i])))
            end
        end
    end
    return dataseries
end

function PlasticityCalibratinator.plotdata_insert!(::Type{Bammann1990Modeling}, modelcalibration)
    # ax, dataseries, Scale_MPa
    Plot_ISVs = []
    for i in range(1, modelcalibration[].modeldata.nsets)
        # println(i)
        scatter!(modelcalibration[].ax,    @lift(Point2f.($(modelcalibration[].dataseries[1][i]).x, $(modelcalibration[].dataseries[1][i]).y)),
            color=i, colormap=:viridis, colorrange=(1, modelcalibration[].modeldata.nsets), label="Data - " * modelcalibration[].modeldata.test_cond["Name"][i])
        lines!(modelcalibration[].ax,      @lift(Point2f.($(modelcalibration[].dataseries[2][i]).x, $(modelcalibration[].dataseries[2][i]).y)),
            color=i, colormap=:viridis, colorrange=(1, modelcalibration[].modeldata.nsets), label="VM Model - " * modelcalibration[].modeldata.test_cond["Name"][i])
        if !isempty(Plot_ISVs)
            scatter!(modelcalibration[].ax_isv,    @lift(Point2f.($(modelcalibration[].dataseries[3][i]).x, $(modelcalibration[].dataseries[3][i]).y)),
                color=i, colormap=:viridis, colorrange=(1, modelcalibration[].modeldata.nsets), label=(s->L"$\alpha$ - %$(s)")(modelcalibration[].modeldata.test_cond["Name"][i]))
            lines!(modelcalibration[].ax_isv,      @lift(Point2f.($(modelcalibration[].dataseries[4][i]).x, $(modelcalibration[].dataseries[4][i]).y)),
                color=i, colormap=:viridis, colorrange=(1, modelcalibration[].modeldata.nsets), label=(s->L"$\kappa$ - %$(s)")(modelcalibration[].modeldata.test_cond["Name"][i]))
            # scatter(modelcalibration[].ax_isv,     @lift(Point2f.($(dataseries[5][i]).x, $(dataseries[5][i]).y)),
            #     color=i, colormap=:viridis , label="\$total\$ - " * bcj.test_cond["Name"][i]))
            # lines(modelcalibration[].ax_isv,       @lift(Point2f.($(dataseries[6][i]).x, $(dataseries[6][i]).y)),
            #     color=i, colormap=:viridis , label="\$S_{11}\$ - " * bcj.test_cond["Name"][i]))
        end
    end
end

# function PlasticityCalibratinator.plotdata_straincontrolkernel(::Type{Bammann1990Modeling}, temp, epsrate, emax, incnum, params)::NamedTuple
#     # println("Setup: emax for set ", i, " = ", emax)
#     # println("θ=", temp, ", ϵ̇=", epsrate, ", ϵₙ=", emax, ", N=", incnum)
#     loading = BCJMetalStrainControl(temp, epsrate, emax, incnum, params)
#     history = kernel(Bammann1990Modeling, loading)[3]
#     # println("Solved: emax for set ", i, " = ", maximum(jc_history.ϵ))
#     return (ϵ=history.ϵ, σ=history.σ, σvM=history.σ)
# end

# @inline function PlasticityCalibratinator.plotdata_updatekernel(::Type{Bammann1990Modeling}, test_data, test_cond, incnum, params, i)::NamedTuple
#     return plotdata_straincontrolkernel(Bammann1990Modeling,
#         test_cond["Temp"][i], test_cond["StrainRate"][i],
#         maximum(test_data["Data_E"][i]), incnum, params)
# end

function PlasticityCalibratinator.plotdata_update!(::Type{Bammann1990Modeling}, modelcalibration)
    # dataseries, incnum, Scale_MPa
    Plot_ISVs = []
    for (key, val) in modelcalibration[].modeldata.params
        modelcalibration[].modeldata.materialproperties[key] = val
    end
    loadstate = modelcalibration[].modeldata.modelinputs.loading_axial ? :tension : :compression
    @sync @distributed for i in range(1, modelcalibration[].modeldata.nsets)
    # for i in range(1, modelcalibration.modeldata.nsets)
        # sol = jccalibration_kernel(jc.test_data, jc.test_cond, jc.incnum, jc.materialproperties, i)
        sol = plotdata_updatekernel(Bammann1990Modeling, modelcalibration[].modeldata.test_data, modelcalibration[].modeldata.test_cond, modelcalibration[].modeldata.incnum, loadstate, modelcalibration[].modeldata.materialproperties, i)
        modelcalibration[].modeldata.test_data["Model_S"][i]    .= sol.σ
        modelcalibration[].modeldata.test_data["Model_VM"][i]   .= sol.σvM
        modelcalibration[].modeldata.test_data["Model_alph"][i] .= sol.α
        modelcalibration[].modeldata.test_data["Model_kap"][i]  .= sol.κ
        modelcalibration[].dataseries[2][i][].y .= modelcalibration[].modeldata.test_data["Model_VM"][i]
        if !isempty(Plot_ISVs)
            modelcalibration[].dataseries[3][i][].y .= modelcalibration[].modeldata.test_data["Model_alph"][i]
            modelcalibration[].dataseries[4][i][].y .= modelcalibration[].modeldata.test_data["Model_kap"][i]
            # modelcalibration[].dataseries[5][i][].y .= modelcalibration[].modeldata.test_data["Model_tot"][i]
            # modelcalibration[].dataseries[6][i][].y .= modelcalibration[].modeldata.test_data["Model_S"][i]
        end
        for ds in modelcalibration[].dataseries[2:end]
            notify(ds[i])
        end
    end
    return nothing
end