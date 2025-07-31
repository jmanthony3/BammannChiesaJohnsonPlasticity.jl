using Dates
using LinearAlgebra
using LUSE_ENGR701_704_NumericalMethods
using Tensors

include("test_dofextraction.jl")

δ(i,j) = i == j ? 1.0 : 0.0 # helper function
Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)

# ::J2Plasticity{}
# J2Plasticity()

# ::MaterialState{}
# MaterialState()

vonMises(x) = (s = dev(x); sqrt(3.0/2.0 * s ⊡ s))
# vonMises(σ) = (s = dev(σ); sign(tr(σ)) * sqrt(3.0/2.0 * s ⊡ s))

# compute_stress_tangent()

function symmetrize_lower!(K)
    for i in axes(K, 1)
        for j in i+1:size(K, 1)
            K[i,j] = K[j,i]
        end
    end
end

# function vectorstrafe(i, n)
#     j = i % (i != n ? n : (n + 1))
#     return [j, j + n, j + 2n]
# end
# [20250123T1359] (JMA3) - (^) close but no cigar | (v) see https://github.com/Ferrite-FEM/Ferrite.jl/discussions/788
vectorstrafe(grid, dh, i) = vertexdofs(dh, nodeid_to_vertexindex(grid, i))
average(array) = sum(array) / length(array)
percentdifference(a, b) = (a == b ? 0.0 : (100 * abs(a - b) / average([a, b])))

"""
(x, y): Actual value

(x̂, ŷ): Predicted value
"""
rmse((x, y), (x̂, ŷ)) = √(length(x) \ sum((ŷ[map(xᵢ->(yᵢ = findfirst(xᵢ .<= x̂); !isnothing(yᵢ) ? yᵢ : findlast(xᵢ .>= x̂)), x)] - y) .^ 2.0))

function derive_serial(x::T, y::T,
        range::UnitRange{Int64}, n::Int64,
        method::Union{String, Symbol},
        δf::Int64=Δf, δt::Float64=Δt)::Float64 where {T<:Vector{Float64}}
    npoint = if isa(method, String)
        (lowercase(method) == "five" ? 5 : (lowercase(method) == "three" ? 3 : nothing))
    elseif isa(method, Symbol)
        (method == :five ? 5 : (method == :three ? 3 : nothing))
    else
        nothing
    end
    return if last(range) - (first(range) - 1) >= npoint
        if first(range) < n < last(range)
            if n - (first(range) - 1) > 2 && last(range) - (n - 1) > 2
                midpoint(x[range], y[range], (δf*δt), n - (first(range) - 1); method=method)
            else
                midpoint(x[range], y[range], (δf*δt), n - (first(range) - 1); method=:three)
            end
        elseif first(range) == n
            endpoint(x[range], y[range], (δf*δt), :begin; method=method)
        elseif n == last(range)
            endpoint(x[range], y[range], (δf*δt), :end; method=method)
        else
            n1derivative(x[n - 1:n + 1], y[n - 1:n + 1], 2)
        end
    else
        if first(range) < n < last(range)
            n1derivative(x[n - 1:n + 1], y[n - 1:n + 1], 2)
        elseif first(range) == n
            n1derivative(x[n:n + 2], y[n:n + 2], 2)
        elseif n == last(range)
            n1derivative(x[n - 2:n], y[n - 2:n], 2)
        else
            Inf
        end
    end
end

timenow()                   = split(@sprintf("%-23s", Dates.now()), "T")[2]
timenow(t::Dates.DateTime)  = split(@sprintf("%-23s", t), "T")[2]

function timeduration(task_0, task_1)
    year            = Dates.year(task_1) - Dates.year(task_0)
    month           = Dates.month(task_1) - Dates.month(task_0)
    day             = Dates.day(task_1) - Dates.day(task_0)
    hour            = Dates.hour(task_1) - Dates.hour(task_0)
    if hour < 0
        day    -= 1;        hour           += 24
    end
    minute          = Dates.minute(task_1) - Dates.minute(task_0)
    if minute < 0
        hour   -= 1;        minute         += 60
    end
    second          = Dates.second(task_1) - Dates.second(task_0)
    if second < 0
        minute -= 1;        second         += 60
    end
    millisecond     = Dates.millisecond(task_1) - Dates.millisecond(task_0)
    if millisecond < 0
        second -= 1;        millisecond    += 1000
    end
    return @sprintf("%04d-%02d-%02dT%02d:%02d:%02d.%03d", year, month, day, hour, minute, second, millisecond)
end