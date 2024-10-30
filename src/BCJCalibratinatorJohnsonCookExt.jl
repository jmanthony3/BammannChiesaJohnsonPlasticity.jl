abstract type JC end

struct JCStrainControl{T<:AbstractFloat} <: JC
    θ       ::T
    ϵ_dot   ::T
    ϵₙ      ::T
    N       ::Integer
    params  ::Dict{String, T}
end

mutable struct JCCurrentConfiguration{T<:AbstractFloat} <: JC
    N       ::Integer
    θ       ::T
    ϵ       ::T
    ϵ_dot   ::T
    Δϵ      ::T
    Tr      ::T
    Tm      ::T
    er0     ::T
    ϵ⁺      ::T
    θ⁺      ::T
    A       ::T
    B       ::T
    n       ::T
    C       ::T
    m       ::T
    σ       ::T
end

mutable struct JCConfigurationHistory{T<:AbstractFloat} <: JC
    σ::Vector{T}
    ϵ::Vector{T}
end

function Base.:+(x::T, y::T) where {T<:JCConfigurationHistory}
    return JCConfigurationHistory{eltype(x.σ)}(
        hcat(x.σ, y.σ),
        hcat(x.ϵ, y.ϵ .+ x.ϵ[:, end])
    )
end

function Base.copyto!(reference::JCCurrentConfiguration, history::JCConfigurationHistory)
    reference.σ = history.σ[:, end]
    return nothing
end

function record!(history::JCConfigurationHistory, i::Integer, current::JCCurrentConfiguration)
    history.σ[i] = current.σ
    history.ϵ[i] = current.ϵ
    return nothing
end

function jcreferenceconfiguration(::Type{JC}, jc::JCStrainControl)::Tuple{JCCurrentConfiguration, JCCurrentConfiguration, JCConfigurationHistory}
    θ       = jc.θ
    ϵ_dot   = jc.ϵ_dot
    ϵₙ      = jc.ϵₙ
    N       = jc.N
    params  = jc.params
    M       = N + 1
    Tr      = params["Tr"]
    Tm      = params["Tm"]
    er0     = params["er0"]
    A       = params["A"]
    B       = params["B"]
    n       = params["n"]
    C       = params["C"]
    m       = params["m"]
    T       = typeof(float(jc.θ))
    ϵ⁺ = ϵ_dot / er0
    θ⁺ = ( θ - Tr ) / ( Tm - Tr )
    Δϵ = ϵₙ/N
    current = JCCurrentConfiguration{T}(N, θ, 0., ϵ_dot, Δϵ,
        Tr, Tm, er0, ϵ⁺, θ⁺, A, B, n, C, m, 0.)
    history = JCConfigurationHistory{T}(
        Vector{T}(undef, M),
        # [range(0., ϵₙ; length=M)...]
        Vector{T}(undef, M)
    )
    record!(history, 1, current)
    return (current, current, history)
end

# JC stress function
function johnsoncookstress(A, B, ϵ, n, C, ϵ⁺, θ⁺, m)::AbstractFloat
    return ( A + B * ϵ^n ) * ( 1. + C * log(ϵ⁺) ) * ( 1. - θ⁺^m )
end

function solve!(jc::JCCurrentConfiguration{<:AbstractFloat},
        history::JCConfigurationHistory)
    A   = jc.A
    B   = jc.B
    n   = jc.n
    C   = jc.C
    ϵ⁺  = jc.ϵ⁺
    θ⁺  = jc.θ⁺
    m   = jc.m
    # calculate the model stress-strain data
    # for (i, ϵ) ∈ zip(range(2, jc.N + 1), history.ϵ[2:end])
    #     history.σ[i] = johnsoncookstress(A, B, ϵ, n, C, ϵ⁺, θ⁺, m)
    # end
    for i ∈ range(2, jc.N + 1)
        jc.ϵ = jc.ϵ + jc.Δϵ
        jc.σ = johnsoncookstress(A, B, jc.ϵ, n, C, ϵ⁺, θ⁺, m)
        record!(history, i, jc)
    end
    return nothing
end