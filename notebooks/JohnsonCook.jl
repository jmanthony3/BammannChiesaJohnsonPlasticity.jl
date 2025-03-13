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
	Pkg.add("ContinuumMechanicsBase")
	Pkg.add("ComponentArrays")
	Pkg.add("StructArrays")
	Pkg.add("DocStringExtensions")
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


	
	using PlutoUI
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

# ╔═╡ 6bdac1c7-2ee0-41da-8240-d595393d6805
begin
	using StructArrays
	using DocStringExtensions
	
	abstract type AbstractJCModel         <: ContinuumMechanicsBase.AbstractMaterialModel end
	abstract type AbstractJCTest{T}       <: ContinuumMechanicsBase.AbstractMaterialTest end
	
	struct JCStrainControl{T1<:Integer, T2<:AbstractFloat} <: ContinuumMechanicsBase.AbstractMaterialTest
	    θ   ::T2    # applied temperature
	    ϵ̇   ::T2    # applied strain rate
	    ϵₙ  ::T2    # final strain
	    N   ::T1    # number of strain increments
	end
	
	struct JCDataEntry{T, S}
	    ϵ::Vector{T}
	    σ::Vector{S}
	end
	
	struct JCUniaxialTest{T, S} <: AbstractJCTest{T}
	    data::StructVector
	    name::String
	    """
	    $(SIGNATURES)
	
	    Creates an object storing results from a uniaxial test of a hyperelatic  material.
	
	    # Arguments:
	    - `ϵ₁`: Vector of experimental, uniaxial strains
	    - `σ₁`: Vector of experimental, uniaxial stresses (optional)
	    - `name`: string for the name of the test
	    - `incompressible`: `true` if the material can be assumed to be incompressible.
	    """
	    function JCUniaxialTest(ϵ₁, σ₁; name, incompressible = true)
	        @assert length(ϵ₁) == length(σ₁) "Inputs must be the same length"
	        if incompressible
	            # ϵ₂ = ϵ₃ = @. sqrt(1 / ϵ₁)
	            ϵ₂ = ϵ₃ = @. -0.499ϵ₁
	        else
	            ϵ₂ = ϵ₃ = Vector{eltype(ϵ₁)}(undef, length(ϵ₁))
	        end
	        ϵ = collect.(zip(ϵ₁, ϵ₂, ϵ₃))
	        σ = collect.(zip(σ₁))
	        data = StructArray{JCDataEntry}((ϵ, σ))
	        new{eltype(eltype(ϵ)), eltype(eltype(σ))}(data, name)
	    end
	    function JCUniaxialTest(ϵ₁; name, incompressible = true)
	        if incompressible
	            # ϵ₂ = ϵ₃ = @. sqrt(1 / ϵ₁)
	            ϵ₂ = ϵ₃ = @. -0.499ϵ₁
	        else
	            ϵ₂ = ϵ₃ = Vector{eltype(ϵ₁)}(undef, length(ϵ₁))
	        end
	        ϵ = collect.(zip(ϵ₁, ϵ₂, ϵ₃))
	        σ = collect.(zip(Vector{eltype(ϵ₁)}(undef, length(ϵ₁))))
	        data = StructArray{JCDataEntry}((ϵ, σ))
	        new{eltype(eltype(ϵ)), eltype(eltype(σ))}(data, name)
	    end
	end
	
	struct JC{T<:AbstractFloat} <: AbstractJCModel
	    θ   ::T         # applied temperature
	    ϵ̇   ::T         # strain rate (effective)
	    ϵₙ  ::T         # final strain
	    N   ::Integer   # number of strain increments
	    Tr  ::T         # reference temperature
	    Tm  ::T         # melting temperature
	    er0 ::T         # reference strain rate
	    ϵ⁺  ::T
	    θ⁺  ::T
	    Δϵ  ::T         # strain increment
	end
	
	function JC(jc::JCStrainControl, Tr::T, Tm::T, er0::T) where {T<:AbstractFloat}
	    ϵ⁺  = jc.ϵ̇ / er0
	    θ⁺  = ( jc.θ - Tr ) / ( Tm - Tr )
	    Δϵ  = jc.ϵₙ / jc.N # strain increment
	    return JC{T}(jc.θ, jc.ϵ̇, jc.ϵₙ, jc.N, Tr, Tm, er0, ϵ⁺, θ⁺, Δϵ)
	end
	
	
	"""
	Function to get a full stress-strain curve (and ISV values)
	
	params = material constants
	
	istate: 1 = tension, 2 = torsion
	
	**no damage in this model**
	"""
	function update(model::JC, ϵ, (;
	            A, B, n, C, m
	        ))
	    return (#=[=#   A   +   (  B * ( ϵ ^ n )  )     #=]=#) * (#=[=#
	        1.0 + ( C * log(model.ϵ⁺) )                 #=]=#) * (#=[=#
	        1.0 - ( model.θ⁺ ^ m )                      #=]=#)
	end
	
	function ContinuumMechanicsBase.predict(
	            ψ   ::JC{T}, # , S},
	            test::AbstractJCTest{T},
	            p;
	            kwargs...,
	        ) where {T<:AbstractFloat} # , S<:SymmetricTensor{2, 3}}
	    M = ψ.N + 1
	    σ     = 0.0 # deviatoric stress
	    ϵₚ    = 0.0 # plastic strain
	    ϵ⃗ = []
	    σ⃗ = []
	    push!(ϵ⃗, ϵₚ)
	    push!(σ⃗, σ)
	    for i ∈ range(2, M)
	        ϵₚ += ψ.Δϵ
	        σ = update(ψ, ϵₚ, p)
	        push!(ϵ⃗, ϵₚ)
	        push!(σ⃗, σ)
	    end
	    return (data=(ϵ=hcat(ϵ⃗...), σ=hcat(σ⃗...)),)
	end
	
	
	function parameters(::M) where {M<:ContinuumMechanicsBase.AbstractMaterialModel} end
	
	function parameter_bounds(::M, ::Any) where {M<:ContinuumMechanicsBase.AbstractMaterialModel}
	    lb = nothing
	    ub = nothing
	    return (lb = lb, ub = ub)
	end
	
	function parameter_bounds(
	            ψ       ::M,
	            tests   ::Vector{Any},
	        ) where {M<:ContinuumMechanicsBase.AbstractMaterialModel}
	    bounds = map(Base.Fix1(parameter_bounds, ψ), tests)
	    lbs = getfield.(bounds, :lb)
	    ubs = getfield.(bounds, :ub)
	    if !(eltype(lbs) <: Nothing)
	        lb_ps = fieldnames(eltype(lbs))
	        lb = map(p -> p .=> maximum(getfield.(lbs, p)), lb_ps) |> NamedTuple
	    else
	        lb = nothing
	    end
	    if !(eltype(ubs) <: Nothing)
	        ub_ps = fieldnames(eltype(ubs))
	        ub = map(p -> p .=> minimum(getfield.(ubs, p)), ub_ps) |> NamedTuple
	    else
	        ub = nothing
	    end
	    return (lb = lb, ub = ub)
	end
	
	parameters(::JC) = (:A, :B, :n, :C, :m)
	
	"""
	$(SIGNATURES)
	
	Creates an `OptimizationProblem` for use in [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) to find the optimal parameters.
	
	# Arguments:
	- `ψ`: material model to use
	- `test` or `tests`: A single or vector of hyperelastics tests to use when fitting the parameters
	- `u₀`: Initial guess for parameters
	- `ps`: Any additional parameters for calling predict
	- `adb`: Select differentiation type from [`ADTypes.jl`](https://github.com/SciML/ADTypes.jl). The type is automatically applied to the type of AD applied to the Optimization Problem also.
	- `loss`: Loss function from [`LossFunctions.jl`](https://github.com/JuliaML/LossFunctions.jl)
	"""
	function JCPlasticityProblem end
	
	function JCPlasticityProblem(
	    ψ   ::JC{T}, # , S},
	    test::JCUniaxialTest{T},
	    u0;
	    ad_type,
	    ui,
	    loss    = L2DistLoss(),
	    lb      = parameter_bounds(ψ, test).lb,
	    ub      = parameter_bounds(ψ, test).ub,
	    int     = nothing,
	    lcons   = nothing,
	    ucons   = nothing,
	    sense   = nothing,
	    kwargs...,
	) where {T<:AbstractFloat} #, S<:SymmetricTensor{2, 3, T}}
	    function f(ps, p)
	        ψ, test, qs, loss, ad_type, kwargs = p
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
	        pred = predict(ψ, test, g(ps, qs); ad_type, kwargs...)
	        resϵ = [first(x) for x in eachcol(pred.data.ϵ)]
	        testϵ = [first(x) for x in test.data.ϵ]
	        s = collect([[x...] for x in eachcol(pred.data.σ)[[findlast(x .>= resϵ) for x in testϵ]]])
	        res = map(i -> loss.(only(i[1]), only(i[2])), zip(s, test.data.σ)) |> mean
	        @show res
	        return res
	    end
	
	    u0 = ComponentVector(u0)
	    if !isnothing(lb) && !isnothing(ub)
	        lb = ComponentVector(lb)
	        ub = ComponentVector(ub)
	    elseif !isnothing(lb)
	        lb = ComponentVector(lb)
	        ub = u0 .* Inf
	    elseif !isnothing(ub)
	        ub = ComponentVector(ub)
	        lb = u0 .* -Inf
	    else
	        ub = u0 .* Inf
	        lb = u0 .* -Inf
	    end
	
	    model_ps = parameters(ψ)
	    for p in model_ps
	        if !isnothing(lb)
	            if (u0[p] < lb[p])
	                @error "Parameter $p = $(u0[p]) is less than lower bound of $(lb[p])"
	                return nothing
	            end
	        end
	        if !isnothing(ub)
	            if (u0[p] > ub[p])
	                @error "Parameter $p = $(u0[p]) is greater than upper bound of $(ub[p])"
	                return nothing
	            end
	        end
	    end
	
	    func = OptimizationFunction(f, ad_type)
	    # Check for Bounds
	    p = (ψ, test, ui, loss, ad_type, kwargs)
	    return OptimizationProblem(func, u0, p; lb, ub, int, lcons, ucons, sense)
	end
end

# ╔═╡ d534bf54-4c83-43d6-a62c-8e4a34f8f74d
md"""
# Johnson-Cook Plasticity Calibration
This notebook can be used to calibrate the constants for the empirical Johnson-Cook plasticity model.
This echos the work of the neighboring notebook for interacting with and calibrating the Bammann-Chiesa-Johnson (BCJ) plasticity model which relies on the the `BammannChiesaJohnsonPlasticity.jl` package; however, this notebook for Johnson-Cook includes a cell that defines all the necessary types/structures, constructors, and functions for method dispatching on the functions defined in `ContinuumMechanicsBase.jl`.
That is, this notebook could be adapted for the rapid development and testing of other material models: e. g. if one wanted to implement their own variation of the Johnson-Cook of BCJ models.
What follows is an example of loading experimental data for _insert test conditions_ and constructing the appropriate BCJ model from test conditions and material properties.

## Initialize Project Environment
"""

# ╔═╡ 741432aa-fbab-4a40-962f-fbd9b5d9166c
md"""
Secondly, we define all the necessary types/structures, constructors, and functions to utilize the method dispatching of `ContinuumMechanicsBase.jl` and perform model calibration.
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
	test = JCUniaxialTest( # (x, y) data from experiment for calibration
		df_Tension_e002_295[!, "Strain"],
		df_Tension_e002_295[!, "Stress"],
		name="exp")
	jc_loading = JCStrainControl( # loading conditions
		295.0, 											# temperature
		2e-3, 											# strain-rate
		float(last(df_Tension_e002_295[!, "Strain"])), 	# final strain
		200) 											# number of strain increments
	nothing
end

# ╔═╡ d6b8bf04-e1fc-41d8-93af-345953f03040
md"""
Now we define some material properties for the desired model.
"""

# ╔═╡ b63e916b-4601-4b61-97ae-9aa07515050c
begin
	Tr = 295.0
	Tm = 1793.0
	er0 = 1.0
	nothing
end

# ╔═╡ a1e992cc-9792-4457-a653-e76d3c47c1da
md"""
Construct the model type given the loading conditions and material properties.
"""

# ╔═╡ 1b83b3e8-9593-483b-a690-fe06aa48aeb5
ψ = JC(jc_loading, Tr, Tm, er0)

# ╔═╡ bd66c9a7-cf0a-4d34-884b-f369722801a8
md"""
Now we can make a group of sliders for the pre-defined model `parameters`.
"""

# ╔═╡ 45ed6284-590e-40ee-93f2-439f264fa032
p0 = ComponentVector(
    A   = 835.014717457143,
    B   = 810.9199208476191,
    n   = 0.5152083333333339,
    C   = 0.004992526238095374,
    m   = 1.1721130952380956
)

# ╔═╡ 2494657a-bdaa-48c5-8209-a36585697975
@bind p parameters_sliders(String.(collect(parameters(ψ))), p0)

# ╔═╡ d4836c95-8b9d-4c0e-bcf3-29abdc551967
p

# ╔═╡ 65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
begin
	res = ContinuumMechanicsBase.predict(ψ, test, p)
	plt = scatter(df_Tension_e002_295[!, "Strain"], df_Tension_e002_295[!, "Stress"], label="exp")
	scatter!(plt, [only(x) for x in eachcol(res.data.ϵ)], [only(x) for x in eachcol(res.data.σ)], label="Johnson-Cook")
end

# ╔═╡ 22a08ebd-2461-4625-8f9b-3ec72cbb5a05
@bind p_checkboxes confirm(MultiCheckBox(String.(collect(parameters(ψ)))))

# ╔═╡ df492d79-2a80-4fb2-ad59-f57f4e2b99e9
# ╠═╡ show_logs = false
begin
	q = parameters_selection(ComponentVector(p), p_checkboxes)
	prob = JCPlasticityProblem(ψ, test, p; ad_type=AutoForwardDiff(), ui=q)
	sol = solve(prob, LBFGS())
	calib = ContinuumMechanicsBase.predict(ψ, test, sol.u)
	scatter!(plt, [only(x) for x in eachcol(calib.data.ϵ)], [only(x) for x in eachcol(calib.data.σ)], label="JC (Calib.)")
end

# ╔═╡ Cell order:
# ╟─d534bf54-4c83-43d6-a62c-8e4a34f8f74d
# ╟─5cacf487-3916-4b7a-8fbf-04c8b4c9a6d9
# ╟─741432aa-fbab-4a40-962f-fbd9b5d9166c
# ╟─6bdac1c7-2ee0-41da-8240-d595393d6805
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
# ╠═65d0598f-fd0b-406b-b53c-3e8b5c4b3d40
# ╠═22a08ebd-2461-4625-8f9b-3ec72cbb5a05
# ╠═df492d79-2a80-4fb2-ad59-f57f4e2b99e9
