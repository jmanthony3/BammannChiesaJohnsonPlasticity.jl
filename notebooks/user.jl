### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
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
	Pkg.add(url="https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl.git", rev="main")
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

# ╔═╡ 5cc1d59a-8722-4bb9-b64b-47a62dfcdeb1
include("user-functions.jl")

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

# ╔═╡ 156a860c-e8a5-4dd8-b234-0a0e4419b5a5
md"""
## Instantiate Calibration Model
### Load Experimental Data
Next, we load the desired `.csv` file and configure the type of material test to calibrate against the experimental data...
"""

# ╔═╡ 398fa1e3-1d11-4285-ad23-b11a4d8628c5
df_Tension_e002_295 = CSV.read("../test/Data_Tension_e0002_T295.csv", DataFrame;
	header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])

# ╔═╡ ba3e98a7-9088-48bf-abeb-110d458b3297
md"""
Using the experimental data, we construct the loading conditions for the initial state for calibration.
"""

# ╔═╡ bd3a90e7-8896-4553-bbd8-bf72c8f60eaf
begin
	test = BCJMetalUniaxialTest( # (x, y) data from experiment for calibration
		df_Tension_e002_295[!, "Strain"],
		df_Tension_e002_295[!, "Stress"] .* 1e6,
		name="exp")
	Ω = BCJMetalStrainControl( # loading conditions
		295.0, 											# temperature
		2e-3, 											# strain-rate
		float(last(df_Tension_e002_295[!, "Strain"])), 	# final strain
		200, 											# number of strain increments
		:tension) 										# load direction
	nothing
end

# ╔═╡ d6b8bf04-e1fc-41d8-93af-345953f03040
md"""
Now we define some material properties for the desired model.
"""

# ╔═╡ b63e916b-4601-4b61-97ae-9aa07515050c
begin
	G 	= 159e9 	# shear modulus [Pa]
	μ 	= 77e9 		# bulk modulus [Pa]
	nothing
end

# ╔═╡ a1e992cc-9792-4457-a653-e76d3c47c1da
md"""
Construct the model type given the loading conditions and material properties.
"""

# ╔═╡ 1b83b3e8-9593-483b-a690-fe06aa48aeb5
ψ = Bammann1993Failure(Ω, μ)

# ╔═╡ bd66c9a7-cf0a-4d34-884b-f369722801a8
md"""
Now we can make a group of sliders for the pre-defined model `parameters`.
"""

# ╔═╡ 45ed6284-590e-40ee-93f2-439f264fa032
p0 = ComponentVector(
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
	m̄ = 1.0,
)

# ╔═╡ 2494657a-bdaa-48c5-8209-a36585697975
@bind p parameters_sliders(String.(collect(parameters(ψ))), p0)

# ╔═╡ d4836c95-8b9d-4c0e-bcf3-29abdc551967
p

# ╔═╡ 65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
begin
	res = ContinuumMechanicsBase.predict(ψ, test, p)
	plt = scatter(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"], label="exp",
		xlabel="True Strain (ϵ) [mm/mm]",
		ylabel="True Stress (σ) [MPa]")
	plot!(plt, [first(x) for x in eachcol(res.data.ϵ)], [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6, label=@sprintf(
			"Bammann1993Failure (RMSE:%.3f)", rmse(
					(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
					([first(x) for x in eachcol(res.data.ϵ)], [vonMises(x) for x in eachcol(res.data.σ)] ./ 1e6))
			),
		linecolor=:blue
	)
end

# ╔═╡ 22a08ebd-2461-4625-8f9b-3ec72cbb5a05
@bind p_checkboxes confirm(MultiCheckBox(String.(collect(parameters(ψ)))))

# ╔═╡ df492d79-2a80-4fb2-ad59-f57f4e2b99e9
begin
	q = parameters_selection(ComponentVector(p), p_checkboxes)
	prob = ContinuumMechanicsBase.MaterialOptimizationProblem(ψ, test, p, parameters(ψ), AutoForwardDiff(), L2DistLoss(), ui=q)
	sol = solve(prob, LBFGS())
	calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
	plot!(deepcopy(plt), [first(x) for x in eachcol(calib.data.ϵ)], [vonMises(x) for x in eachcol(calib.data.σ)] ./ 1e6, label=@sprintf(
			"Bammann1993Failure (RMSE:%.3f, K:%d, T:%.3f [s])", rmse(
				(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"]),
				([first(x) for x in eachcol(calib.data.ϵ)], [vonMises(x) for x in eachcol(calib.data.σ)] ./ 1e6)),
			sol.stats.iterations, sol.stats.time),
		linecolor=:blue,
		linestyle=:dash)
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
# ╠═5cc1d59a-8722-4bb9-b64b-47a62dfcdeb1
# ╟─156a860c-e8a5-4dd8-b234-0a0e4419b5a5
# ╟─398fa1e3-1d11-4285-ad23-b11a4d8628c5
# ╟─ba3e98a7-9088-48bf-abeb-110d458b3297
# ╠═bd3a90e7-8896-4553-bbd8-bf72c8f60eaf
# ╟─d6b8bf04-e1fc-41d8-93af-345953f03040
# ╠═b63e916b-4601-4b61-97ae-9aa07515050c
# ╟─a1e992cc-9792-4457-a653-e76d3c47c1da
# ╠═1b83b3e8-9593-483b-a690-fe06aa48aeb5
# ╟─bd66c9a7-cf0a-4d34-884b-f369722801a8
# ╠═45ed6284-590e-40ee-93f2-439f264fa032
# ╠═2494657a-bdaa-48c5-8209-a36585697975
# ╠═d4836c95-8b9d-4c0e-bcf3-29abdc551967
# ╟─65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
# ╠═22a08ebd-2461-4625-8f9b-3ec72cbb5a05
# ╠═df492d79-2a80-4fb2-ad59-f57f4e2b99e9
# ╠═ac027691-ae47-4450-b9d6-b814b5be79d5
