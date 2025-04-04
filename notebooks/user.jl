### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 5cacf487-3916-4b7a-8fbf-04c8b4c9a6d9
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".") # activate project in current directory
	Pkg.add(url="https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl.git", rev="cho2019unified")
	Pkg.add("ContinuumMechanicsBase")
	Pkg.add("ComponentArrays")
	Pkg.add("CSV")
	Pkg.add("DataFrames")
	Pkg.add("DocStringExtensions")
	Pkg.add("FiniteDiff")
	Pkg.add("ForwardDiff")
	Pkg.add("Optimization")
	Pkg.add("OptimizationOptimJL")
	Pkg.add("LossFunctions")
	Pkg.add("Plots")
	Pkg.add("Pluto")
	Pkg.add("PlutoUI")
	Pkg.add("Printf")



	using PlutoUI
	import PlutoUI: combine

	using BammannChiesaJohnsonPlasticity
	using ContinuumMechanicsBase
	using ComponentArrays
	using CSV, DataFrames
	using DocStringExtensions
	using FiniteDiff
	import ForwardDiff
	using Optimization, OptimizationOptimJL, LossFunctions
	using Plots
	using Printf

	# include("user-functions.jl")



	function parameters_sliders(parameters::Vector, values::ComponentVector)
		return PlutoUI.combine() do Child
			
			inputs = [
				md""" $(parameter): $(
					Child(parameter, Slider(value .* logrange(1e-3, 1e3, length=1001), default=value))
				)"""
				
				for (parameter, value) in zip(parameters, values)
			]
			
			md"""
			#### Sliders for Coarse Adjustment of Material Constants
			$(inputs)
			"""
		end
	end

	function parameters_checkboxes(parameters::Vector)
		return PlutoUI.combine() do Child
			
			inputs = [
				md"""$(parameter): $(CheckBox())"""
				
				for parameter in parameters
			]
			
			md"""
			#### Checkboxes for Fine Adjustment (Optimization) of Selected Material Constants
			$(inputs)
			"""
		end
	end

	function parameters_selection(parameters::ComponentVector, checkboxes::Vector)
		if !isempty(checkboxes)
			for key in checkboxes
				parameters[key] = NaN
			end
		end
		return ComponentVector(parameters)
	end
end

# ╔═╡ d534bf54-4c83-43d6-a62c-8e4a34f8f74d
md"""
# Bammann-Chiesa-Johnson Plasticity Calibration
This notebook can be used to calibrate the constants for Internal State Variables (ISVs) in the Julian implementation of the Bammann-Chiesa-Johnson (BCJ) plasticity model.
The `BammannChiesaJohnsonPlasticity.jl` package was modeled after the `Hyperelastics.jl` package implementing the `ContinuumMechanicsBase.jl` package for common function signatures.
Therefore, the `BCJPlasticity.jl` package is fully capable of performing point-simulator predictions for a variety of uniaxial loading conditions at given temperatures and strain rates.
With the `Optimization.jl` package from SciML, calibrations for BCJ model constants may also be performed.
This notebook expands on the `BCJPlasticity.jl` package with sliders, checkboxes, and other widgets from the `PlutoUI.jl` package which adds a layer of interaction with the BCJ plasticity model of choice.
What follows is an example of loading experimental data from a tension test of 4340 stainless steel at room temperature ($295 [K]$) and $2 \times 10^{-3} [mm/mm/s]$ strain rate and constructing the appropriate BCJ model from test conditions and material properties.

## Initialize Project Environment
First, we start by loading the required packages and defining some helper functions.
"""

# ╔═╡ 8e1775ee-a0a0-4538-a606-b758cb6f302c
function ContinuumMechanicsBase.MaterialOptimizationProblem(
    ψs   ::Vector{<:AbstractBCJModel},  # , S},
    tests::Vector{<:AbstractBCJTest},
    u₀,
    model_ps,
    ad_type,
    loss;
    ui,
    lb      = parameter_bounds(first(ψs), first(tests)).lb,
    ub      = parameter_bounds(first(ψs), first(tests)).ub,
    int     = nothing,
    lcons   = nothing,
    ucons   = nothing,
    sense   = nothing,
    kwargs...,
) # where {T<:AbstractFloat} #, S<:SymmetricTensor{2, 3, T}}
    get_data(d, f; columnate=false, concatenate=false) = (z = [[first(x) for x in (columnate ? eachcol(y.data[f]) : y.data[f])] for y in d]; concatenate ? vcat(z...) : z)
    get_idx(d1, d2) = [[findlast(first(x) .>= get_data([d2[i]], :ϵ; columnate=true)...) for x in z] for (i, z) in enumerate([[first(x) for x in y.data.ϵ] for y in values(d1)])]
    function f(ps, p)
        ψs, tests, qs, loss, ad_type, kwargs = p
        function g(ps, qs)
            if !isnothing(qs) && any(!isnan, qs)
                for (name, value) in zip(keys(qs), qs)
                    if !isnan(value)
                        # @show value
                        ps[name] = value
                    end
                end
            end
            return ComponentVector(ps)
        end
        preds = [ContinuumMechanicsBase.predict(ψ, test, g(ps, qs); ad_type, kwargs...) for (ψ, test) in zip(ψs, tests)]
        # @show preds
        # resϵ = [first(x) for x in eachcol(pred.data.ϵ)]
        # testϵ = [first(x) for x in test.data.ϵ]
        # ϵ = get_data(tests, ϵ)
        # ϵ̂ = get_data(preds, ϵ; col=true)
        # # resϵ = [x[1, 1] for x in pred.data.ϵ]
        # # testϵ = [x[1, 1] for x in test.data.ϵ]
        # s = collect([[x...] for x in eachcol(pred.data.σ)[[findlast(x .>= resϵ) for x in testϵ]]])
        # # s = collect([[x...] for x in pred.data.σ[[findlast(x .>= resϵ) for x in testϵ]]])
        s = vcat([x[y] for (x, y) in zip(get_data(preds, :ϵ; columnate=true), get_idx(tests, preds))]...)
        ŝ = vcat([[first(x) for x in y.data.ϵ] for y in values(tests)]...)
        # @show length(s) == length(ŝ)
        res = map(i -> loss.(only(i[1]), only(i[2])), zip(s, ŝ)) |> mean
        # @show res # uncomment for testing
        return res
    end

    u₀ = ComponentVector(u₀)
    # pb = ContinuumMechanicsBase.parameter_bounds(ψ, test)
    # lb, ub = pb.lb, pb.ub
    if !isnothing(lb) && !isnothing(ub)
        lb = ComponentVector(lb)
        ub = ComponentVector(ub)
    elseif !isnothing(lb)
        lb = ComponentVector(lb)
        ub = u₀ .* Inf
    elseif !isnothing(ub)
        ub = ComponentVector(ub)
        lb = u₀ .* -Inf
    else
        ub = u₀ .* Inf
        lb = u₀ .* -Inf
    end

    # model_ps = ContinuumMechanicsBase.parameters(ψ)
    for p in model_ps
        if !isnothing(lb)
            if (u₀[p] < lb[p])
                @error "Parameter $p = $(u₀[p]) is less than lower bound of $(lb[p])"
                return nothing
            end
        end
        if !isnothing(ub)
            if (u₀[p] > ub[p])
                @error "Parameter $p = $(u₀[p]) is greater than upper bound of $(ub[p])"
                return nothing
            end
        end
    end

    func = OptimizationFunction(f, ad_type)
    # Check for Bounds
    p = (ψs, tests, ui, loss, ad_type, kwargs)
    return OptimizationProblem(func, u₀, p; lb, ub, int, lcons, ucons, sense)
end

# ╔═╡ 5cc1d59a-8722-4bb9-b64b-47a62dfcdeb1
# ╠═╡ disabled = true
#=╠═╡
include("user-functions.jl")
  ╠═╡ =#

# ╔═╡ 156a860c-e8a5-4dd8-b234-0a0e4419b5a5
md"""
## Instantiate Calibration Model
### Load Experimental Data
Next, we load the desired `.csv` file and configure the type of material test to calibrate against the experimental data...
"""

# ╔═╡ 398fa1e3-1d11-4285-ad23-b11a4d8628c5
# df_Tension_e002_295 = CSV.read("../test/Data_Tension_e0002_T295.csv", DataFrame;
# 	header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
df_Fig4a = CSV.read("Cho2019UnifiedStaticDynamic-Fig4a.csv", DataFrame;
        header=true, skipto=3, delim=',', types=Float64)

# ╔═╡ ba3e98a7-9088-48bf-abeb-110d458b3297
md"""
Using the experimental data, we construct the loading conditions for the initial state for calibration.
"""

# ╔═╡ d6b8bf04-e1fc-41d8-93af-345953f03040
md"""
Now we define some material properties for the desired model.
"""

# ╔═╡ b63e916b-4601-4b61-97ae-9aa07515050c
begin
	ϵ̇ = 4e-4
	K = 159e9   # bulk modulus [Pa]
	μ = 77e9    # shear modulus [Pa]
	nothing
end

# ╔═╡ a1e992cc-9792-4457-a653-e76d3c47c1da
md"""
Construct the model type given the loading conditions and material properties.
"""

# ╔═╡ 1b83b3e8-9593-483b-a690-fe06aa48aeb5
begin
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
	    # @show (4(i - 1) + 1, 4(i - 1) + 2), θ_str, ϵ̇, last(x), 4length(x)
	    tests[θ_str] = BCJMetalUniaxialTest(x, y, name="$(θ_flt)K")
	    domains[θ_str] = BCJMetalStrainControl(θ_flt, ϵ̇, last(x), 4length(x), :tension)
	    models[θ_str] = Bammann1990Modeling(domains[θ_str], μ)
	end
	tests = sort(tests; rev=false)
	domains = sort(domains; rev=false)
	models = sort(models; rev=false)
end

# ╔═╡ bd66c9a7-cf0a-4d34-884b-f369722801a8
md"""
Now we can make a group of sliders for the pre-defined model `parameters`.
"""

# ╔═╡ 45ed6284-590e-40ee-93f2-439f264fa032
p0 = ComponentVector(
    C₁ = 9.1402e10,
    C₂ = 258.417,
    C₃ = 1.62805e8,
    C₄ = 363.053,
    C₅ = 1.28544,
    C₆ = 236.047,
    C₇ = 1.04959e-6,
    C₈ = 0.0920373,
    C₉ = 4.07014e-10,
    C₁₀ = 1000.0,
    C₁₁ = 7.07701e-12,
    C₁₂ = 18.6325,
    C₁₃ = 5.07815e-12,
    C₁₄ = 38.7783,
    C₁₅ = 3.77314e7,
    C₁₆ = 0.0111427,
    C₁₇ = 7.87311e6,
    C₁₈ = 0.0155747,
)

# ╔═╡ 2494657a-bdaa-48c5-8209-a36585697975
@bind p parameters_sliders(String.(collect(parameters(first(models)[2]))), p0)

# ╔═╡ d4836c95-8b9d-4c0e-bcf3-29abdc551967
p

# ╔═╡ 65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
begin
    plt = plot(xlims=(0, 1), ylims=(0, Inf), widen=1.06)
	for (i, (θ, ψ)) in enumerate(models)
        test = tests[θ]
        res = ContinuumMechanicsBase.predict(ψ, test, p)
        # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
        scatter!(plt, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ],
                markercolor=i,
                label="$(θ)K:Exp",
            )
        plot!(plt, [first(x) for x in eachcol(res.data.ϵ)], [vonMises(x) for x in eachcol(res.data.σ)],
                linecolor=i,
                label="$(θ)K:Model",
            )
    end
    plt
end

# ╔═╡ 22a08ebd-2461-4625-8f9b-3ec72cbb5a05
@bind p_checkboxes confirm(MultiCheckBox(String.(collect(parameters(first(models)[2])))))

# ╔═╡ df492d79-2a80-4fb2-ad59-f57f4e2b99e9
begin
	pltq = plot(xlims=(0, 1), ylims=(0, Inf), widen=1.06)
	q = parameters_selection(ComponentVector(p), p_checkboxes)
	# prob = ContinuumMechanicsBase.MaterialOptimizationProblem(ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(), ui=q)
	# sol = solve(prob, LBFGS())
	# calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
	# plot!(deepcopy(plt), [first(x) for x in eachcol(calib.data.ϵ)], [vonMises(x) for x in eachcol(calib.data.σ)] ./ 1e6, label=@sprintf(
	# 		"Bammann1993Failure (RMSE:%.3f, K:%d, T:%.3f [s])", rmse(
	# 			(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
	# 			([first(x) for x in eachcol(calib.data.ϵ)], [vonMises(x) for x in eachcol(calib.data.σ)] ./ 1e6)),
	# 		sol.stats.iterations, sol.stats.time),
	# 	linecolor=:blue,
	# 	linestyle=:dash)
	prob = ContinuumMechanicsBase.MaterialOptimizationProblem(
	    collect(Bammann1990Modeling, values(models)),
	    collect(BCJMetalUniaxialTest, values(tests)),
	    p,
	    parameters(first(values(models))),
	    AutoForwardDiff(),
	    L2DistLoss();
	    ui=q)
	sol = solve(prob, LBFGS())
	for (i, (θ, ψ)) in enumerate(models)
        test = tests[θ]
        calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
        # @show [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6
		scatter!(pltq, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ],
                markercolor=i,
                label="$(θ)K:Exp",
            )
        plot!(pltq, [first(x) for x in eachcol(calib.data.ϵ)], [vonMises(x) for x in eachcol(calib.data.σ)],
                linecolor=i,
                label="$(θ)K:Calib",
            )
    end
	pltq
end

# ╔═╡ ac027691-ae47-4450-b9d6-b814b5be79d5
@show sol.retcode; i, r = 1, deepcopy(q); for (key, value) in zip(keys(p), q)
	if isnan(value)
		r[key] = sol.u[i]
		@printf("\t%s = %.9f,\n", key, r[key])
	end
	global i += 1
end; r

# ╔═╡ Cell order:
# ╟─d534bf54-4c83-43d6-a62c-8e4a34f8f74d
# ╠═5cacf487-3916-4b7a-8fbf-04c8b4c9a6d9
# ╠═8e1775ee-a0a0-4538-a606-b758cb6f302c
# ╠═5cc1d59a-8722-4bb9-b64b-47a62dfcdeb1
# ╟─156a860c-e8a5-4dd8-b234-0a0e4419b5a5
# ╟─398fa1e3-1d11-4285-ad23-b11a4d8628c5
# ╟─ba3e98a7-9088-48bf-abeb-110d458b3297
# ╟─d6b8bf04-e1fc-41d8-93af-345953f03040
# ╠═b63e916b-4601-4b61-97ae-9aa07515050c
# ╟─a1e992cc-9792-4457-a653-e76d3c47c1da
# ╠═1b83b3e8-9593-483b-a690-fe06aa48aeb5
# ╟─bd66c9a7-cf0a-4d34-884b-f369722801a8
# ╠═45ed6284-590e-40ee-93f2-439f264fa032
# ╠═2494657a-bdaa-48c5-8209-a36585697975
# ╠═d4836c95-8b9d-4c0e-bcf3-29abdc551967
# ╠═65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
# ╠═22a08ebd-2461-4625-8f9b-3ec72cbb5a05
# ╠═df492d79-2a80-4fb2-ad59-f57f4e2b99e9
# ╠═ac027691-ae47-4450-b9d6-b814b5be79d5
