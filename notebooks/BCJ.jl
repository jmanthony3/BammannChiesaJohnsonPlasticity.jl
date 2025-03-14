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
	# this package is only needed until BCJPlasticity.jl#CMB is merged/registered
	Pkg.add(url="https://github.com/jmanthony3/PlasticityBase.jl.git")
	Pkg.add(url="https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl.git", rev="CMB")
	Pkg.add("ContinuumMechanicsBase")
	Pkg.add("ComponentArrays")
	Pkg.add("CSV")
	Pkg.add("DataFrames")
	Pkg.add("FiniteDiff")
	Pkg.add("ForwardDiff")
	Pkg.add("Optimization")
	Pkg.add("OptimizationOptimJL")
	Pkg.add("LossFunctions")
	Pkg.add("Plots")
	Pkg.add("Pluto")
	Pkg.add("PlutoUI")


	
	using BammannChiesaJohnsonPlasticity, PlutoUI
	import PlutoUI: combine

	using ContinuumMechanicsBase
	using ComponentArrays
	using CSV, DataFrames
	using FiniteDiff
	import ForwardDiff
	using Optimization, OptimizationOptimJL, LossFunctions
	using Plots



	function parameters_sliders(parameters::Vector, values::ComponentVector)
		return PlutoUI.combine() do Child
			
			inputs = [
				md""" $(parameter): $(
					Child(parameter, Slider(value .* logrange(0.1, 10.0, length=1000), default=value))
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

	md"""
	First, we start by loading the required packages and defining some helper functions.
	"""
end

# ╔═╡ d534bf54-4c83-43d6-a62c-8e4a34f8f74d
md"""
# Bammann-Chiesa-Johnson Plasticity Calibration
This notebook can be used to calibrate the constants for Internal State Variables (ISVs) in the Julian implementation of the Bammann-Chiesa-Johnson (BCJ) plasticity model.
The `BammannChiesaJohnsonPlasticity.jl` package was modelled after the `Hyperelastics.jl` package implementing the `ContinuumMechanicsBase.jl` package for standard function signatures.
Therefore, the `BCJPlasticity.jl` package is fully capable of performing point-simulator calibrations for BCJ model constants via the `Optimization.jl` package from SciML.
What follows is an example of loading experimental data for _insert test conditions_ and constructing the appropriate BCJ model from test conditions and material properties.

## Initialize Project Environment
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
	bcj_loading = BCJMetalStrainControl( # loading conditions
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
ψ = DK(bcj_loading, μ)

# ╔═╡ bd66c9a7-cf0a-4d34-884b-f369722801a8
md"""
Now we can make a group of sliders for the pre-defined model `parameters`.
"""

# ╔═╡ 45ed6284-590e-40ee-93f2-439f264fa032
p0 = ComponentVector(
	C₁ 	= 35016459.896579415, 		C₂ 	= 323.93342698083165, 	# V
	C₃ 	= 500340419.8337271, 		C₄ 	= 143.08381901004486, 	# Y
	C₅ 	= 4.101775377562497, 		C₆ 	= 271.0245526, 			# f
	C₇ 	= 1.0834796217232945e-06, 	C₈ 	= 1023.6278003945317, 	# r_s
	C₉ 	= 2358205093.844017, 		C₁₀ = 676421.9935474312, 	# h
	C₁₁ = 1.3465080192134937e-10, 	C₁₂ = 98.35671405000001, 	# r_d
	C₁₃ = 2.533629073577668e-09, 	C₁₄ = 403.2291451343492, 	# R_s
	C₁₅ = 1159915808.5023918, 		C₁₆ = 959557.0948847248, 	# H
	C₁₇ = 6.204370386543724e-12, 	C₁₈ = 203.95288011132806, 	# R_s
	C₁₉ = 1e-10, 					C₂₀ = 1e-10 				# Y_adj
)

# ╔═╡ 2494657a-bdaa-48c5-8209-a36585697975
@bind p parameters_sliders(String.(collect(parameters(ψ))), p0)

# ╔═╡ d4836c95-8b9d-4c0e-bcf3-29abdc551967
p

# ╔═╡ 65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
begin
	res = ContinuumMechanicsBase.predict(ψ, test, p)
	plt = scatter(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"] .* 1e6, label="exp")
	scatter!(plt, [first(x) for x in eachcol(res.data.ϵ)], [symmetricvonMises(x) for x in eachcol(res.data.σ)], label="DK")
end

# ╔═╡ 22a08ebd-2461-4625-8f9b-3ec72cbb5a05
@bind p_checkboxes confirm(MultiCheckBox(String.(collect(parameters(ψ)))))

# ╔═╡ df492d79-2a80-4fb2-ad59-f57f4e2b99e9
# ╠═╡ show_logs = false
begin
	q = parameters_selection(ComponentVector(p), p_checkboxes)
	prob = ContinuumMechanicsBase.MaterialOptimizationProblem(ψ, test, p, AutoForwardDiff(), L2DistLoss(), ui=q)
	sol = solve(prob, LBFGS())
	calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
	scatter!(plt, [first(x) for x in eachcol(calib.data.ϵ)], [symmetricvonMises(x) for x in eachcol(calib.data.σ)], label="DK (Calib.)")
end

# ╔═╡ Cell order:
# ╟─d534bf54-4c83-43d6-a62c-8e4a34f8f74d
# ╟─5cacf487-3916-4b7a-8fbf-04c8b4c9a6d9
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
