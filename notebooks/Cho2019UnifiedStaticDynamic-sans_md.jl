begin
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

	# ╔═╡ 5624f9cd-d46d-4dfd-ad46-e6185c120eba
	# ╠═╡ show_logs = false
	# begin
	# 	using Pkg
	# 	Pkg.activate(".") # activate project in current directory
	# 	Pkg.add(url="https://github.com/jmanthony3/BammannChiesaJohnsonPlasticity.jl.git", rev="cho2019unified")
	# 	Pkg.add("ContinuumMechanicsBase")
	# 	Pkg.add("ComponentArrays")
	# 	Pkg.add("CSV")
	# 	Pkg.add("DataFrames")
	# 	Pkg.add("DataInterpolations")
	# 	Pkg.add("DocStringExtensions")
	# 	Pkg.add("FiniteDiff")
	# 	Pkg.add("ForwardDiff")
	# 	Pkg.add("Optimization")
	# 	Pkg.add("OptimizationOptimJL")
	# 	Pkg.add("LossFunctions")
	# 	Pkg.add("Plots")
	# 	Pkg.add("Pluto")
	# 	Pkg.add("PlutoUI")
	# 	Pkg.add("Printf")
	# 	Pkg.add("HypertextLiteral")
	# 	Pkg.add("Distributed")



	# 	using PlutoUI
	# 	import PlutoUI: combine

	# 	using BammannChiesaJohnsonPlasticity
	# 	using ContinuumMechanicsBase
	# 	using ComponentArrays
	# 	using CSV, DataFrames
	# 	using DocStringExtensions
	# 	using FiniteDiff
	# 	import ForwardDiff
	# 	using DataInterpolations
	# 	using Optimization, OptimizationOptimJL, LossFunctions
	# 	using Plots
	# 	using Printf
	# 	using HypertextLiteral
	# 	using Distributed
	# 	begin
	# 		function parameters_sliders(parameter_groups, value_groups::Vector)
	# 			return PlutoUI.combine() do Child
	# 				inputs = []
	# 				j, k = 1, 0
	# 				for i in range(1, length(value_groups))
	# 					for value_group in value_groups[i]
	# 						k += length(value_group)
	# 						parameter_group = parameter_groups[j:k]
	# 						# input_group = [md" $(parameter): $(
	# 						# 			   Child(parameter, Slider(value .* logrange(1e-3, 1e3, length=1001), default=value, show_value=true))
	# 						# 			   ) " for (parameter, value) in zip(parameter_group, value_group)]
	# 						# Scrubbable(0.80 : 0.01 : 1.00, format=".0%", prefix="you are 🌝 ", suffix=" cool")
	# 						input_group = [md""" $(parameter): $(
	# 									Child(parameter, Scrubbable(value .* logrange(1e-3, 1e3; length=1001); default=value, format=".03e")))
	# 									) """ for (parameter, value) in zip(parameter_group, value_group)]
	# 						push!(inputs, md" $(input_group...) ")
	# 						j += length(value_group)
	# 					end
	# 				end
					
	# 				# md"""
	# 				# #### Sliders for Coarse Adjustment of Material Constants
	# 				# $(inputs)
	# 				# """
	# 				md" $inputs "
	# 			end
	# 		end

	# 		function parameters_checkboxes(parameters::Vector)
	# 			return PlutoUI.combine() do Child
					
	# 				inputs = [
	# 					md"""$(parameter): $(CheckBox())"""
						
	# 					for parameter in parameters
	# 				]
					
	# 				md"""
	# 				#### Checkboxes for Fine Adjustment (Optimization) of Selected Material Constants
	# 				$(inputs)
	# 				"""
	# 			end
	# 		end

	# 		function parameters_selection(parameters::ComponentVector, checkboxes::Vector)
	# 			if !isempty(checkboxes)
	# 				for key in checkboxes
	# 					parameters[key] = NaN
	# 				end
	# 			end
	# 			return ComponentVector(parameters)
	# 		end

	# 		# # https://fonsp-disorganised-mess.netlify.app/layout
	# 		# function aside(x)
	# 		# 	@htl("""
	# 		# 		<style>
					
					
	# 		# 		@media (min-width: calc(700px + 30px + 300px)) {
	# 		# 			aside.plutoui-aside-wrapper {
	# 		# 				position: absolute;
	# 		# 				right: -11px;
	# 		# 				width: 0px;
	# 		# 			}
	# 		# 			aside.plutoui-aside-wrapper > div {
	# 		# 				width: 300px;
	# 		# 			}
	# 		# 		}
	# 		# 		</style>
					
	# 		# 		<aside class="plutoui-aside-wrapper">
	# 		# 		<div>
	# 		# 		$(x)
	# 		# 		</div>
	# 		# 		</aside>
					
	# 		# 		""")
	# 		# end
	# 		# aside(embed_display(p))

	# 		# Pkg.add("PrettyTables")
	# 		# using PrettyTables
	# 		# df = DataFrame(rand(100, 5), :auto) # Example DataFrame

	# 		# html_table = pretty_table(String, df, backend = Val(:html))
			
	# 		# scrollable_table = @htl("""
	# 		# <div style="max-height: 400px; overflow-y: auto; border: 1px solid #ccc;">
	# 		#     $(HTML(html_table))
	# 		# </div>
	# 		# """)
			
	# 		# scrollable_table
	# 	end
	# end

	default(fontfamily="Computer Modern", linewidth=1, framestyle=:box)

	# ╔═╡ 5cc1d59a-8722-4bb9-b64b-47a62dfcdeb1
	# include("Cho2019UnifiedStaticDynamic-functions+new.jl")
	include("Cho2019UnifiedStaticDynamic-functions.jl")
end
begin
	# ╔═╡ 398fa1e3-1d11-4285-ad23-b11a4d8628c5
	df_Fig4a = CSV.read("Cho2019UnifiedStaticDynamic-Fig4a.csv", DataFrame;
		header=true, delim=',', skipto=3, types=Float64)
	# df_Fig4a = CSV.read("data-SS.csv", DataFrame;
	# 	header=true, delim=',', skipto=3, types=Float64)

	# ╔═╡ b63e916b-4601-4b61-97ae-9aa07515050c
	begin
		n 	= 2.0
		ω₀ 	= 3.6e4
		R 	= 8.31446261815324 # universal gas constant
		# R 	= 8.314 # universal gas constant
		E⁺ 	= 82.0e3
		# E⁺ 	= 82.0
		V⁺ 	= 0.0
		z 	= 0.65
		# z 	= 0.0
		# d₀ 	= 10.0 # μm (Ghauri et al., 1990)
		d₀ 	= 62.0 # μm (Tanner et al., 1990)
		η₀ 	= 0.0
		Kic = 1000.0
		# Kic = 50.0
		𝒹 	= 0.0
		# 𝒹 	= 2.0e-5
		𝒻 	= 0.001
		R₀ 	= 0.0
		# R₀ 	= 1e-6
		nothing
	end

	# ╔═╡ bd3a90e7-8896-4553-bbd8-bf72c8f60eaf
	# ╠═╡ show_logs = false
	begin
		ϵ̇, N = 4e-4, 1000
		tests = Dict()
		domains = Dict()
		models = Dict()
		for (i, θ) in enumerate((298, 407, 475, 509, 542, 559, 576, 610, 678, 814))
			# if θ == 298
				θ_str = match(r"(.*)K(.*)", names(df_Fig4a)[4(i - 1) + 1])[1]
				θ_flt = parse(Float64, θ_str)
				x = filter(!ismissing, df_Fig4a[!, 4(i - 1) + 1])
				idx_sort = sortperm(x)
				x = x[idx_sort]
				y = filter(!ismissing, df_Fig4a[!, 4(i - 1) + 2])[idx_sort] .* 1e6
				@show (4(i - 1) + 1, 4(i - 1) + 2), θ_str, ϵ̇, ceil(1.01last(x); sigdigits=2), 4length(x)
				tests[θ_str] = BCJMetalUniaxialTest(x, y, name="$(θ_flt)K")
				# domains[θ_str] = BCJMetalStrainControl(θ_flt, ϵ̇, last(x), 4length(x), :compression)
				domains[θ_str] = BCJMetalStrainControl(θ_flt, ϵ̇, ceil(1.01last(x); sigdigits=2), N, :compression)
				models[θ_str] = Cho2019UnifiedStaticDynamic(domains[θ_str], n, ω₀, E⁺, V⁺, R, d₀, z, Kic, 𝒹, 𝒻, η₀, R₀)
			# end
		end
		
		tests = sort(tests; rev=false)
		domains = sort(domains; rev=false)
		models = sort(models; rev=false)
		nothing
	end

	# ╔═╡ 45ed6284-590e-40ee-93f2-439f264fa032
	p0 = [
		[
			ComponentVector(C₁ = 5.637, C₂ = 112.6,),
			ComponentVector(C₃ = 8.378, C₄ = 324.9,),
			ComponentVector(C₅ = 2.971, C₆ = 2548.0,),
			# ComponentVector(C₁ = 5.5, C₂ = 112.6,),
			# ComponentVector(C₃ = 5.278, C₄ = 324.9,),
			# ComponentVector(C₅ = 1970, C₆ = 2548.0,),
		],
		[
			ComponentVector(Pₖ₁ = 0.0, Pₖ₂ = 0.0, Pₖ₃ = 0.0,),
		],
		[
			ComponentVector(C₇ = 0.1345, C₈ = 351.1, C₂₁ = 0.0,),
			ComponentVector(C₉ = 0.02869, C₁₀ = 0.0, C₂₂ = 0.0,),
			ComponentVector(C₁₁ = 0.02928, C₁₂ = 4337.0, C₂₃ = 0.0,),
			# ComponentVector(C₇ = 0.1845, C₈ = 351.1, C₂₁ = 0.0,),
			# ComponentVector(C₉ = 0.0255, C₁₀ = 0.0, C₂₂ = 0.0,),
			# ComponentVector(C₁₁ = 40.28, C₁₂ = 4337.0, C₂₃ = 0.0,),
		],
		[
			ComponentVector(C₁₃ = 0.05098, C₁₄ = 476.6, C₂₄ = 0.0,),
			ComponentVector(C₁₅ = 0.006924, C₁₆ = 0.0, C₂₅ = 0.0,),
			ComponentVector(C₁₇ = 2.487, C₁₈ = 7611.0, C₂₆ = 0.0,),
			ComponentVector(NK = 2.0,),
			# ComponentVector(C₁₃ = 0.05098, C₁₄ = 476.6, C₂₄ = 0.0,),
			# ComponentVector(C₁₅ = 7.524e-3, C₁₆ = 0.0, C₂₅ = 0.0,),
			# ComponentVector(C₁₇ = 54.87, C₁₈ = 7611.0, C₂₆ = 0.0,),
			# ComponentVector(NK = 2.0,),
			# ComponentVector(C₁₃ = 0.05098, C₁₄ = 476.6, C₂₄ = 0.0,),
			# ComponentVector(C₁₅ = 7.524e-3, C₁₆ = 0.0, C₂₅ = 0.0,),
			# ComponentVector(C₁₇ = 54.87, C₁₈ = 7611.0, C₂₆ = 0.0,),
			# ComponentVector(NK = 2.0,),
		],
		[
			ComponentVector(C₁₉ = 0.0, C₂₀ = 0.0,),
		],
		[
			ComponentVector(Cx1 = 1.78e6, Cx2 = 7.806e3, Cdp = 0.0,),
			ComponentVector(Cx3 = 5.401e4, Cx4 = 8943.0, Csp = 0.0,),
			ComponentVector(Cx5 = 5.0, Cxa = 0.8052, Cxb = 3.68, Cxc = 4.485,),
		],
		[
			ComponentVector(Cg1 = 7.41e4, Cg2 = 0.8826, Cg3 = 1.185e-3,),
			# ComponentVector(Cg1 = 7.41e4, Cg2 = 0.8826, Cg3 = 1.185e-3, Cg4 = 0.0,),
			# ComponentVector(Cg1 = 7.4e4, Cg2 = 0.8800, Cg3 = 1.2e-3,),
			# ComponentVector(Cg1 = 7.4e4, Cg2 = 0.8800, Cg3 = 1.2e-2, Cg4 = 0.0,),
		],
		[
			ComponentVector(Dc = 0.99, a = 0.0, b = 0.0, c = 0.0),
			# ComponentVector(a = 0.0, b = 0.0, c = 3.3e4,),
		],
		[
			ComponentVector(Cnuc = 0.0, Tnuc = 0.0, nn = 0.3, Tgrw = 0.0,),
			# ComponentVector(Cnuc = 1.0e15, Tnuc = 10.0, nn = 0.3, Tgrw = 0.0,),
		],
		[
			ComponentVector(kr1 = 0.0, krt = 0.0, kr2 = 0.0,
			kr3 = 0.0, kp1 = 0.0, kpt = 0.0, kp2 = 0.0,),
			# ComponentVector(kr1 = 7e-32, krt = 5e3, kr2 = 3.5,
			# kr3 = 4.1e2, kp1 = 1.4e-27, kpt = 2.5e3, kp2 = 2.8,),
		]
	]

	# ╔═╡ 53926f5c-e18c-4cb6-b062-bb965ec41769
	slider_ui = @bind p parameters_sliders(String.(parameters(first(models)[2])), p0);

	plt = plot(xlims=(0, 1), ylims=(0, Inf), legendposition=:best, widen=1.06)
	plt_α = deepcopy(plt)
	plt_d = deepcopy(plt)
	plt_X = deepcopy(plt)
	# @sync @distributed for (i, (θ, ψ)) in collect(enumerate(models))
	for (i, (θ, ψ)) in collect(enumerate(models))
		test = tests[θ]
		prediction = ContinuumMechanicsBase.predict(ψ, test, p; imat=1, iYS=0, iREXmethod=6, iGSmethod=5)
		# prediction = ContinuumMechanicsBase.predict(ψ, test, p; imat=1, iYS=0, iREXmethod=6, iGSmethod=3)
		@info i, θ
		@show vonMises(last(eachcol(prediction.data.σ)))
		@show last(prediction.data.α), last(prediction.data.κ)
		@show last(prediction.data.ϕ), last(prediction.data.X), last(prediction.data.d)
		begin
			# # i, plt = 1, plot(xlims=(0, 1), ylims=(0, Inf), legendposition=:outerright, widen=1.06)
			# i, plt = 1, plot(xlims=(0, 1), ylims=(0, Inf), legendposition=:best, widen=1.06)
			# plt_α = deepcopy(plt)
			# plt_d = deepcopy(plt)
			# plt_X = deepcopy(plt)
			scatter!(plt, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ] ./ 1e6,
					markercolor=i,
					label="$(θ)K:Exp",
				)
			plot!(plt, [first(x) for x in eachcol(prediction.data.ϵ)], [vonMises(x) for x in eachcol(prediction.data.σ)],
					linecolor=i,
					label="$(θ)K:Model",
				)
			ϵ = [first(x) for x in test.data.ϵ]
			# @show ϵ
			σ = [first(x) for x in test.data.σ] ./ 1e6
			# @show σ
			ϵ̂ = [first(x) for x in eachcol(prediction.data.ϵ)]
			# @show ϵ̂
			σ̂ = [vonMises(x) for x in eachcol(prediction.data.σ)]
			# @show prediction.data.σ
			# @show [first(x) for x in eachcol(prediction.data.σ)]
			# @show σ̂
			# s = linear_interpolation(ϵ̂, σ̂, extrapolation_bc=Line()).(ϵ)
			s = CubicSpline(σ̂, ϵ̂; extrapolation=ExtrapolationType.Linear).(ϵ)
			println(rmse((ϵ, σ), (ϵ, s)))
			plot!(plt_α, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.α),
					linecolor=i,
					linestyle=:dash,
					label="$(θ)K:Model,α",
				)
			plot!(plt_α, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.κ),
					linecolor=i,
					label="$(θ)K:Model,κ",
				)
			# ? [20250723T1114] (JMA3): For some reason, `d` immediately goes to 0 which doesn't seem right.
			plot!(plt_d, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.d),
					linecolor=i,
					label="$(θ)K:Model,d",
				)
			plot!(plt_X, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.X),
					linecolor=i,
					label="$(θ)K:Model,X",
				)
			plot!(plt_X, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.ϕ),
					linecolor=i,
					linestyle=:dash,
					label="$(θ)K:Model,ϕ",
				)
		end
	end
	plot(plt, plt_α, plt_d, plt_X; layout=(2, 2), legend=false)
	# plot(plt, legend=false)
end
plot(plt, legendposition=:outerright)
# plot(plt, plt_α, plt_d, plt_X; layout=(2, 2), legend=false)
# # prediction = ContinuumMechanicsBase.predict(ψ, test, p; imat=1, iYS=2, iREXmethod=3, iGSmethod=1)
# # prediction = ContinuumMechanicsBase.predict(ψ, test, p; imat=1, iYS=0, iREXmethod=3, iGSmethod=4)
# prediction = ContinuumMechanicsBase.predict(ψ, test, p; imat=1, iYS=0, iREXmethod=3, iGSmethod=5)
# # prediction = ContinuumMechanicsBase.predict(ψ, test, p; imat=1, iYS=0, iREXmethod=3, iGSmethod=3)

# begin
# 	# i, plt = 1, plot(xlims=(0, 1), ylims=(0, Inf), legendposition=:outerright, widen=1.06)
# 	i, plt = 1, plot(xlims=(0, 1), ylims=(0, Inf), legendposition=:best, widen=1.06)
# 	plt_α = deepcopy(plt)
# 	plt_d = deepcopy(plt)
# 	plt_X = deepcopy(plt)
# 	scatter!(plt, [first(x) for x in test.data.ϵ], [first(x) for x in test.data.σ] ./ 1e6,
# 			markercolor=i,
# 			label="$(θ)K:Exp",
# 		)
# 	plot!(plt, [first(x) for x in eachcol(prediction.data.ϵ)], [vonMises(x) for x in eachcol(prediction.data.σ)],
# 			linecolor=i,
# 			label="$(θ)K:Model",
# 		)
# 	ϵ = [first(x) for x in test.data.ϵ]
# 	# @show ϵ
# 	σ = [first(x) for x in test.data.σ] ./ 1e6
# 	# @show σ
# 	ϵ̂ = [first(x) for x in eachcol(prediction.data.ϵ)]
# 	# @show ϵ̂
# 	σ̂ = [vonMises(x) for x in eachcol(prediction.data.σ)]
# 	# @show prediction.data.σ
# 	# @show [first(x) for x in eachcol(prediction.data.σ)]
# 	# @show σ̂
# 	# s = linear_interpolation(ϵ̂, σ̂, extrapolation_bc=Line()).(ϵ)
# 	s = CubicSpline(σ̂, ϵ̂; extrapolation=ExtrapolationType.Linear).(ϵ)
# 	println(rmse((ϵ, σ), (ϵ, s)))
# 	plot!(plt_α, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.α),
# 			linecolor=i,
# 			linestyle=:dash,
# 			label="$(θ)K:Model,α",
# 		)
# 	plot!(plt_α, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.κ),
# 			linecolor=i,
# 			label="$(θ)K:Model,κ",
# 		)
# 	# ? [20250723T1114] (JMA3): For some reason, `d` immediately goes to 0 which doesn't seem right.
# 	plot!(plt_d, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.d),
# 			linecolor=i,
# 			label="$(θ)K:Model,d",
# 		)
# 	plot!(plt_X, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.X),
# 			linecolor=i,
# 			label="$(θ)K:Model,X",
# 		)
# 	plot!(plt_X, [first(x) for x in eachcol(prediction.data.ϵ)], vec(prediction.data.ϕ),
# 			linecolor=i,
# 			linestyle=:dash,
# 			label="$(θ)K:Model,ϕ",
# 		)
# 	plot(plt, plt_α, plt_d, plt_X; layout=(2, 2))
# end