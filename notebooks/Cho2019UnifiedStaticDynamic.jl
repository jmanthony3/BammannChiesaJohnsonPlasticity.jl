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

# ╔═╡ 5cc1d59a-8722-4bb9-b64b-47a62dfcdeb1
include("Cho2019UnifiedStaticDynamic-functions.jl")

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
# df_Tension_e002_295 = CSV.read("../test/Data_Tension_e0002_T295.csv", DataFrame;
# 	header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
df_Fig4a = CSV.read("Cho2019UnifiedStaticDynamic-Fig4a.csv", DataFrame;
	header=true, delim=',', skipto=3, types=Float64)

# ╔═╡ d6b8bf04-e1fc-41d8-93af-345953f03040
md"""
Now we define some material properties for the desired model.
"""

# ╔═╡ b63e916b-4601-4b61-97ae-9aa07515050c
begin
	n 	= 2.0
	ω₀ 	= 3.6e4
	R 	= 8.31446261815324 # universal gas constant
	E⁺ 	= 82.0e3
	z 	= 0.65
	d₀ 	= 10.0 # μm (Ghauri et al., 1990)
	η₀ 	= 0.0
	Kic = 1000.0
	𝒹 	= 0.0
	𝒻 	= 0.001
	R₀ 	= 0.0
	nothing
end

# ╔═╡ ba3e98a7-9088-48bf-abeb-110d458b3297
md"""
Using the experimental data, we construct the loading conditions for the initial state for calibration.
"""

# ╔═╡ a1e992cc-9792-4457-a653-e76d3c47c1da
md"""
Construct the model type given the loading conditions and material properties.
"""

# ╔═╡ bd3a90e7-8896-4553-bbd8-bf72c8f60eaf
begin
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
	    models[θ_str] = Cho2019Unified(domains[θ_str], E⁺, E⁺, R, d₀, Kic, 𝒹, 𝒻, η₀, R₀)
	end
	
	tests = sort(tests; rev=false)
	domains = sort(domains; rev=false)
	models = sort(models; rev=false)
	nothing
end

# ╔═╡ bd66c9a7-cf0a-4d34-884b-f369722801a8
md"""
Now we can make a group of sliders for the pre-defined model `parameters`.
"""

# ╔═╡ 45ed6284-590e-40ee-93f2-439f264fa032
p0 = ComponentVector(
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
	n = n,
	ω₀ = ω₀,
	Cg1 = 7.41e4,
	Cg2 = 0.8826,
	Cg3 = 1.185e-3,
	z = z,
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
	    collect(Cho2019Unified, values(models)),
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
# ╠═5cc1d59a-8722-4bb9-b64b-47a62dfcdeb1
# ╟─156a860c-e8a5-4dd8-b234-0a0e4419b5a5
# ╟─398fa1e3-1d11-4285-ad23-b11a4d8628c5
# ╟─d6b8bf04-e1fc-41d8-93af-345953f03040
# ╠═b63e916b-4601-4b61-97ae-9aa07515050c
# ╟─ba3e98a7-9088-48bf-abeb-110d458b3297
# ╟─a1e992cc-9792-4457-a653-e76d3c47c1da
# ╠═bd3a90e7-8896-4553-bbd8-bf72c8f60eaf
# ╟─bd66c9a7-cf0a-4d34-884b-f369722801a8
# ╠═45ed6284-590e-40ee-93f2-439f264fa032
# ╠═2494657a-bdaa-48c5-8209-a36585697975
# ╠═d4836c95-8b9d-4c0e-bcf3-29abdc551967
# ╟─65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
# ╠═22a08ebd-2461-4625-8f9b-3ec72cbb5a05
# ╠═df492d79-2a80-4fb2-ad59-f57f4e2b99e9
# ╠═ac027691-ae47-4450-b9d6-b814b5be79d5
