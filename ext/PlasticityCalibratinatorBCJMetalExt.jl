function PlasticityCalibratinator.plotdata_straincontrolkernel(T::Type{<:BCJMetal}, temp, epsrate, emax, incnum, loadtype, params)::NamedTuple
    kS          = 1     # default tension component
    if loadtype == :torsion
        kS      = 4     # select torsion component
    end
    bcj_loading = BCJMetalStrainControl(temp, epsrate, emax, incnum, loadtype, params)
    bcj_history = kernel(T, bcj_loading)[3]
    ϵ__         = bcj_history.ϵ__
    σ__         = bcj_history.σ__
    α__         = bcj_history.α__
    κ           = bcj_history.κ

    # pull only the relevant (tension/torsion) strain being evaluated:
    ϵ       = ϵ__[kS, :]
    σ       = σ__[kS, :]
    σvM     = symmetricvonMises(σ__)
    α       = α__[kS, :]
    return (ϵ=ϵ, σ=σ, σvM=σvM, α=α, κ=κ)
end

@inline function PlasticityCalibratinator.plotdata_updatekernel(T::Type{<:BCJMetal}, test_data, test_cond, incnum, loadtype, params, i)::NamedTuple
    return plotdata_straincontrolkernel(T, test_cond["Temp"][i], test_cond["StrainRate"][i],
        maximum(test_data["Data_E"][i]), incnum, loadtype, params)
end