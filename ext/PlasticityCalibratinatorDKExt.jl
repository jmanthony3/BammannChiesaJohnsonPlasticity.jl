using CSV
using DataFrames
using Distributed
using GLMakie
using LaTeXStrings

set_theme!(theme_latexfonts())

materialprops_DK = Dict(
    # Comment,For Calibration with vumat
    "C01"       => 35016459.896579415,
    "C02"       => 323.93342698083165,
    "C03"       => 500340419.8337271,
    "C04"       => 143.08381901004486,
    "C05"       => 4.101775377562497,
    "C06"       => 271.0245526,
    "C07"       => 1.0834796217232945e-06,
    "C08"       => 1023.6278003945317,
    "C09"       => 2358205093.844017,
    "C10"       => 676421.9935474312,
    "C11"       => 1.3465080192134937e-10,
    "C12"       => 98.35671405000001,
    "C13"       => 2.533629073577668e-09,
    "C14"       => 403.2291451343492,
    "C15"       => 1159915808.5023918,
    "C16"       => 959557.0948847248,
    "C17"       => 6.204370386543724e-12,
    "C18"       => 203.95288011132806,
    "C19"       => 1e-10,
    "C20"       => 1e-10,
    "Bulk Mod"  => 159000000000.0,
    "Shear Mod" => 77000000000.0,
)
materialconsts_DK = collect(
    "C01"       => 35016459.896579415,
    "C02"       => 323.93342698083165,
    "C03"       => 500340419.8337271,
    "C04"       => 143.08381901004486,
    "C05"       => 4.101775377562497,
    "C06"       => 271.0245526,
    "C07"       => 1.0834796217232945e-06,
    "C08"       => 1023.6278003945317,
    "C09"       => 2358205093.844017,
    "C10"       => 676421.9935474312,
    "C11"       => 1.3465080192134937e-10,
    "C12"       => 98.35671405000001,
    "C13"       => 2.533629073577668e-09,
    "C14"       => 403.2291451343492,
    "C15"       => 1159915808.5023918,
    "C16"       => 959557.0948847248,
    "C17"       => 6.204370386543724e-12,
    "C18"       => 203.95288011132806,
    "C19"       => 1e-10,
    "C20"       => 1e-10,
)

PlasticityCalibratinator.materialproperties(::Type{DK}) = materialprops_DK
PlasticityCalibratinator.materialconstants(::Type{DK}) = materialconsts_DK

PlasticityCalibratinator.characteristicequations(::Type{DK}) = [
    # plasticstrainrate
    L"\dot{\epsilon}_{p} = f(\theta)\sinh\left[ \frac{ \{|\mathbf{\xi}| - \kappa - Y(\theta) \} }{ (1 - \phi)V(\theta) } \right]\frac{\mathbf{\xi}'}{|\mathbf{\xi}'|}\text{, let }\mathbf{\xi}' = \mathbf{\sigma}' - \mathbf{\alpha}'",
    # kinematic hardening
    L"\dot{\mathbf{\alpha}} = h(\theta)\dot{\epsilon}_{p} - [r_{d}(\theta)|\dot{\epsilon}_{p}| + r_{s}(\theta)]|\mathbf{\alpha}|\mathbf{\alpha}",
    # isotropic hardening
    L"\dot{\kappa} = H(\theta)\dot{\epsilon}_{p} - [R_{d}(\theta)|\dot{\epsilon}_{p}| + R_{s}(\theta)]\kappa^{2}",
    # flow rule
    L"F = |\sigma - \alpha| - \kappa - \beta(|\dot{\epsilon}_{p}|, \theta)",
    # initial yield stress beta
    L"\beta(\dot{\epsilon}_{p}, \theta) = Y(\theta) + V(\theta)\sinh^{-1}\left(\frac{|\dot{\epsilon}_{p}|}{f(\theta)}\right)"
]
PlasticityCalibratinator.dependenceequations(::Type{DK})     = [
    L"V = C_{ 1} \mathrm{exp}(-C_{ 2} / \theta)",       # V
    L"Y = C_{ 3} \mathrm{exp}( C_{ 4} / \theta)",       # Y
    L"f = C_{ 5} \mathrm{exp}(-C_{ 6} / \theta)",       # f
    L"r_{d} = C_{ 7} \mathrm{exp}(-C_{ 8} / \theta)",   # rd
    L"h = C_{ 9} - C_{10}\theta",                       # h
    L"r_{s} = C_{11} \mathrm{exp}(-C_{12} / \theta)",   # rs
    L"R_{d} = C_{13} \mathrm{exp}(-C_{14} / \theta)",   # Rd
    L"H = C_{15} - C_{16}\theta",                       # H
    L"R_{s} = C_{17} \mathrm{exp}(-C_{18} / \theta)",   # Rs
    L"Y_{adj}",                                         # Yadj
]
PlasticityCalibratinator.dependencesliders(::Type{DK})       = [
    # V
    [
        (label=L"C_{ 1}", range=range(0., 5materialconsts_DK["C01"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C01"]),
        (label=L"C_{ 2}", range=range(0., 5materialconsts_DK["C02"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C02"]),
    ],
    # Y
    [
        (label=L"C_{ 3}", range=range(0., 5materialconsts_DK["C03"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C03"]),
        (label=L"C_{ 4}", range=range(0., 5materialconsts_DK["C04"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C04"]),
    ],
    # f
    [
        (label=L"C_{ 5}", range=range(0., 5materialconsts_DK["C05"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C05"]),
        (label=L"C_{ 6}", range=range(0., 5materialconsts_DK["C06"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C06"]),
    ],
    # rd
    [
        (label=L"C_{ 7}", range=range(0., 5materialconsts_DK["C07"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C07"]),
        (label=L"C_{ 8}", range=range(0., 5materialconsts_DK["C08"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C08"]),
    ],
    # h
    [
        (label=L"C_{ 9}", range=range(0., 5materialconsts_DK["C09"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C09"]),
        (label=L"C_{10}", range=range(0., 5materialconsts_DK["C10"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C10"]),
    ],
    # rs
    [
        (label=L"C_{11}", range=range(0., 5materialconsts_DK["C11"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C11"]),
        (label=L"C_{12}", range=range(0., 5materialconsts_DK["C12"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C12"]),
    ],
    # Rd
    [
        (label=L"C_{13}", range=range(0., 5materialconsts_DK["C13"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C13"]),
        (label=L"C_{14}", range=range(0., 5materialconsts_DK["C14"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C14"]),
    ],
    # H
    [
        (label=L"C_{15}", range=range(0., 5materialconsts_DK["C15"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C15"]),
        (label=L"C_{16}", range=range(0., 5materialconsts_DK["C16"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C16"]),
    ],
    # Rs
    [
        (label=L"C_{17}", range=range(0., 5materialconsts_DK["C17"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C17"]),
        (label=L"C_{18}", range=range(0., 5materialconsts_DK["C18"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C18"]),
    ],
    # Yadj
    [
        (label=L"C_{19}", range=range(0., 5materialconsts_DK["C19"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C19"]),
        (label=L"C_{20}", range=range(0., 5materialconsts_DK["C10"]; length=1_000), format="{:.3e}", startvalue=materialconsts_DK["C10"]),
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
function PlasticityCalibratinator.modeldata(::Type{DK}, modelinputs::ModelInputs, params)::ModelData
    # files, incnum, params, Scale_MPa
    files = modelinputs.expdatasets
    incnum = modelinputs.incnum
    Scale_MPa = modelinputs.stressscale
    # istate = Int(modelinputs.loading_axial)
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
            sol = plotdata_updatekernel(DK, test_data, test_cond, incnum, loadstate, params, i)
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
    return ModelData{DK}(DK, modelinputs, nsets, test_data, test_cond, params, deepcopy(materialconstants(DK)), deepcopy(materialconstants(DK)), incnum, Scale_MPa)
end

function PlasticityCalibratinator.plotdata_initialize(::Type{DK}, nsets, test_data)
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

function PlasticityCalibratinator.plotdata_insert!(::Type{DK}, modelcalibration)
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

# function PlasticityCalibratinator.plotdata_straincontrolkernel(::Type{DK}, temp, epsrate, emax, incnum, params)::NamedTuple
#     # println("Setup: emax for set ", i, " = ", emax)
#     # println("θ=", temp, ", ϵ̇=", epsrate, ", ϵₙ=", emax, ", N=", incnum)
#     loading = BCJMetalStrainControl(temp, epsrate, emax, incnum, params)
#     history = kernel(DK, loading)[3]
#     # println("Solved: emax for set ", i, " = ", maximum(jc_history.ϵ))
#     return (ϵ=history.ϵ, σ=history.σ, σvM=history.σ)
# end

# @inline function PlasticityCalibratinator.plotdata_updatekernel(::Type{DK}, test_data, test_cond, incnum, params, i)::NamedTuple
#     return plotdata_straincontrolkernel(DK,
#         test_cond["Temp"][i], test_cond["StrainRate"][i],
#         maximum(test_data["Data_E"][i]), incnum, params)
# end

function PlasticityCalibratinator.plotdata_update!(::Type{DK}, modelcalibration)
    Plot_ISVs = []
    # dataseries, incnum, Scale_MPa
    for (key, val) in modelcalibration[].modeldata.params
        modelcalibration[].modeldata.materialproperties[key] = val
    end
    loadstate = modelcalibration[].modeldata.modelinputs.loading_axial ? :tension : :compression
    @sync @distributed for i in range(1, modelcalibration[].modeldata.nsets)
    # for i in range(1, modelcalibration.modeldata.nsets)
        # sol = jccalibration_kernel(jc.test_data, jc.test_cond, jc.incnum, jc.materialproperties, i)
        sol = plotdata_updatekernel(DK, modelcalibration[].modeldata.test_data, modelcalibration[].modeldata.test_cond, modelcalibration[].modeldata.incnum, loadstate, modelcalibration[].modeldata.materialproperties, i)
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

# dora

doramat_DK = collect(
    "Bulk Mod"  => 159000000000.0,
    "Shear Mod" => 77000000000.0,
)
PlasticityCalibratinator.materialdora(::Type{DK}) = doramat_DK
dorarange_DK_max = collect([
    "Bulk Mod"  => 10e9,
    "Shear Mod" => 10e6,
])
PlasticityCalibratinator.doraequations(::Type{DK})     = [L"G", L"\mu"]
PlasticityCalibratinator.dorasliders(::Type{DK})       = [
    (label=L"G",    range=range(doramat_DK["Bulk Mod"] - dorarange_DK_max["Bulk Mod"],      doramat_DK["Bulk Mod"] + dorarange_DK_max["Bulk Mod"];      length=1_000), format="{:.3e}", startvalue=doramat_DK["Bulk Mod"]),
    (label=L"\mu",  range=range(doramat_DK["Shear Mod"] - dorarange_DK_max["Shear Mod"],    doramat_DK["Shear Mod"] + dorarange_DK_max["Shear Mod"];    length=1_000), format="{:.3e}", startvalue=doramat_DK["Shear Mod"]),
]

function PlasticityCalibratinator.modeldora(::Type{DK}, model_inputs::ModelInputs, params, temperatures, strainrates, finalstrains)::ModelData
    # files, incnum, params, Scale_MPa
    nsets = mapreduce(length, *, (temperatures, strainrates, finalstrains))
    incnum = model_inputs.incnum
    Scale_MPa = model_inputs.stressscale
    # istate = Int(modelinputs.loading_axial)
    loadstate = model_inputs.loading_axial ? :tension : :compression
    test_cond   = Dict(
        "StrainRate"    => [],
        "StrainFinal"   => [],
        "Temp"          => [],
        "Name"          => []
    )
    test_data   = Dict(
        "Data_E"        => [],
        "Model_E"       => [],
        "Model_S"       => [],
        "Model_alph"    => [],
        "Model_kap"     => [],
        "Model_VM"      => []
    )
    i = 1
    for temp in temperatures
        for strainrate in strainrates
            for finalstrain in finalstrains
                # store the stress-strain data
                name = "$finalstrain @ ($temp, $strainrate)"
                push!(test_cond["StrainRate"],  (i, strainrate))
                push!(test_cond["StrainFinal"], (i, finalstrain))
                push!(test_cond["Temp"],        (i, temp))
                push!(test_cond["Name"],        (i, name))
                push!(test_data["Data_E"],      (i, range(0, finalstrain; length=incnum+1)))
                println((i, strainrate, finalstrain, temp, name))
                i += 1
            end
        end
    end
    for key in ("StrainRate", "StrainFinal", "Temp", "Name")
        for (i, x) in enumerate(sort(test_cond[key], by=x->first(x)))
            test_cond[key][i] = last(x)
        end
    end
    for key in ("Data_E",)
        for (i, x) in enumerate(sort(test_data[key], by=x->first(x)))
            test_data[key][i] = last(x)
        end
    end
    @sync @distributed for i in range(1, nsets)
    # for i in range(1, nsets)
        # test_data, test_cond, incnum, istate, params, i, ISV_Model
        sol = plotdata_straincontrolkernel(DK, test_cond["Temp"][i], test_cond["StrainRate"][i], test_cond["StrainFinal"][i], incnum, loadstate, params)
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
    return ModelData{DK}(DK, model_inputs, nsets, test_data, test_cond, params, deepcopy(materialconstants(DK)), deepcopy(materialconstants(DK)), incnum, Scale_MPa)
end

function PlasticityCalibratinator.plotdora_initialize(::Type{DK}, nsets, test_data)
    plot_isvs = []
    dataseries = if !isempty(plot_isvs)
        [[], [], []]
    else
        [[]]
    end
    if nsets > 0
        for i in range(1, nsets)
            push!(dataseries[1], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_VM"][i])))
            if !isempty(plot_isvs)
                push!(dataseries[2], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_alph"][i])))
                push!(dataseries[3], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_kap"][i])))
                # push!(dataseries[4], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_tot"][i])))
                # push!(dataseries[5], Observable(DataFrame(x=test_data["Model_E"][i], y=test_data["Model_S"][i])))
            end
        end
    end
    return dataseries
end

function PlasticityCalibratinator.plotdora_insert!(::Type{DK}, model_calibration)
    # ax, dataseries, Scale_MPa
    Plot_ISVs = []
    for i in range(1, model_calibration[].modeldata.nsets)
        # println(i)
        lines!(model_calibration[].ax,      @lift(Point2f.($(model_calibration[].dataseries[1][i]).x, $(model_calibration[].dataseries[1][i]).y)),
            color=i, colormap=:viridis, colorrange=(0, model_calibration[].modeldata.nsets), label="VM Model - " * model_calibration[].modeldata.test_cond["Name"][i])
        if !isempty(Plot_ISVs)
            scatter!(model_calibration[].ax_isv,    @lift(Point2f.($(model_calibration[].dataseries[2][i]).x, $(model_calibration[].dataseries[2][i]).y)),
                color=i, colormap=:viridis, colorrange=(0, model_calibration[].modeldata.nsets), label=(s->L"$\alpha$ - %$(s)")(model_calibration[].modeldata.test_cond["Name"][i]))
            lines!(model_calibration[].ax_isv,      @lift(Point2f.($(model_calibration[].dataseries[3][i]).x, $(model_calibration[].dataseries[3][i]).y)),
                color=i, colormap=:viridis, colorrange=(0, model_calibration[].modeldata.nsets), label=(s->L"$\kappa$ - %$(s)")(model_calibration[].modeldata.test_cond["Name"][i]))
            # scatter(modelcalibration[].ax_isv,     @lift(Point2f.($(dataseries[5][i]).x, $(dataseries[5][i]).y)),
            #     color=i, colormap=:viridis , label="\$total\$ - " * bcj.test_cond["Name"][i]))
            # lines(modelcalibration[].ax_isv,       @lift(Point2f.($(dataseries[6][i]).x, $(dataseries[6][i]).y)),
            #     color=i, colormap=:viridis , label="\$S_{11}\$ - " * bcj.test_cond["Name"][i]))
        end
    end
end

function PlasticityCalibratinator.plotdora_update!(::Type{DK}, model_calibration)
    Plot_ISVs = []
    # dataseries, incnum, Scale_MPa
    for (key, val) in model_calibration[].modeldata.params
        model_calibration[].modeldata.materialproperties[key] = val
    end
    loadstate = model_calibration[].modeldata.modelinputs.loading_axial ? :tension : :compression
    @sync @distributed for i in range(1, model_calibration[].modeldata.nsets)
    # for i in range(1, modelcalibration.modeldata.nsets)
        # sol = jccalibration_kernel(jc.test_data, jc.test_cond, jc.incnum, jc.materialproperties, i)
        # sol = plotdata_updatekernel(DK, model_calibration[].modeldata.test_data, model_calibration[].modeldata.test_cond, model_calibration[].modeldata.incnum, loadstate, model_calibration[].modeldata.materialproperties, i)
        sol = plotdata_straincontrolkernel(DK, model_calibration[].modeldata.test_cond["Temp"][i], model_calibration[].modeldata.test_cond["StrainRate"][i], model_calibration[].modeldata.test_cond["StrainFinal"][i], model_calibration[].modeldata.incnum, loadstate, model_calibration[].modeldata.materialproperties)
        model_calibration[].modeldata.test_data["Model_S"][i]    .= sol.σ
        model_calibration[].modeldata.test_data["Model_VM"][i]   .= sol.σvM
        model_calibration[].modeldata.test_data["Model_alph"][i] .= sol.α
        model_calibration[].modeldata.test_data["Model_kap"][i]  .= sol.κ
        model_calibration[].dataseries[1][i][].y .= model_calibration[].modeldata.test_data["Model_VM"][i]
        if !isempty(Plot_ISVs)
            model_calibration[].dataseries[2][i][].y .= model_calibration[].modeldata.test_data["Model_alph"][i]
            model_calibration[].dataseries[3][i][].y .= model_calibration[].modeldata.test_data["Model_kap"][i]
            # modelcalibration[].dataseries[4][i][].y .= modelcalibration[].modeldata.test_data["Model_tot"][i]
            # modelcalibration[].dataseries[5][i][].y .= modelcalibration[].modeldata.test_data["Model_S"][i]
        end
        for ds in model_calibration[].dataseries
            notify(ds[i])
        end
    end
    return nothing
end