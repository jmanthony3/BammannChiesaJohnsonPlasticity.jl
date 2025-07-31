# CDFEM.jl
# 
# see the following for relevant discussion and example code
#   - von Mises plasticity example (plain program): https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/plasticity/
#   - apply nodal force for unstruct mesh: https://github.com/Ferrite-FEM/Ferrite.jl/discussions/1128
#   - The proper way of setting concentrated loads on the right-hand side vector (node to DoF index): https://github.com/Ferrite-FEM/Ferrite.jl/discussions/788
#   - Confusion of Units and Ability to Apply Fixed Displacement: https://github.com/Ferrite-FEM/Ferrite.jl/discussions/1137#discussioncomment-12019789
# [20250217T1650] (JMA3) - Meeting in Dr. Cho's Office
#   - Ensure that ALL material model constants are same across experimental, calibrated, and FEA implementations
#   - Ensure that strain rate between experimental, calibrated, and FEA data is consistent between timesteps
#   - TODO: create plots of strain rates

using CSV
using DataFrames
using Ferrite
using FerriteGmsh
using LinearAlgebra
# using OhMyThreads
using Plots
using Printf
using SparseArrays
# using TaskLocalValues
using Tensors
using WriteVTK

include("src/functions/preamble.jl")

struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}}
    G   ::T # shear modulus
    K   ::T # bulk modulus
    σ₀  ::T # initial yield limit
    H   ::T # hardening modulus
    Dᵉ  ::S # elastic stiffness tensor
end;

function J2Plasticity(E, ν, σ₀, H)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
    return J2Plasticity(G, K, σ₀, H, Dᵉ)
end;

struct JCPlasticity{T, S <: SymmetricTensor{4, 3, T}}
    G   ::T # shear modulus
    K   ::T # bulk modulus
    σ₀  ::T # initial yield limit
    Dᵉ  ::S # elastic stiffness tensor
end;

function JCPlasticity(E, ν, σ₀)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
    return JCPlasticity(G, K, σ₀, Dᵉ)
end;

struct MaterialState{T, S <: SecondOrderTensor{3, T}}
    # store "converged" values
    ϵ   ::S
    ϵᵖ  ::S # plastic strain
    σ   ::S # stress
    k   ::T # hardening variable
end

function MaterialState()
    return MaterialState(
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0)
end

struct JohnsonCookState{T, S <: SecondOrderTensor{3, T}}
    # store "converged" values
    ϵ   ::S
    ϵᵖ  ::S # plastic strain
    σ   ::S # stress
    k   ::T # hardening variable
end

function JohnsonCookState()
    return JohnsonCookState(
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0)
end

# function compute_stress_tangent(ϵ::SymmetricTensor{2, 3}, material::J2Plasticity, state::MaterialState)
function compute_stress_tangent(Δϵ::SymmetricTensor{2, 3}, material::J2Plasticity, state::MaterialState)
    # unpack some material parameters
    G = material.G
    H = material.H

    # we use (•)ᵗ to denote *trial*-values
    σᵗ = material.Dᵉ ⊡ (ϵ - state.ϵᵖ) # trial-stress
    Δϵ = ϵ - state.ϵ

    # # σᵗ = state.σ + material.Dᵉ ⊡ (ϵ - state.ϵ) # trial-stress
    # # Δϵ = ϵ - state.ϵ
    # σᵗ = state.σ + material.Dᵉ ⊡ Δϵ # trial-stress
    # # σᵗ = material.Dᵉ ⊡ (state.ϵ + Δϵ) # trial-stress

    # sᵗ = dev(σᵗ)         # deviatoric part of trial-stress
    # J₂ = 0.5 * sᵗ ⊡ sᵗ   # second invariant of sᵗ
    # σᵗₑ = sqrt(3.0*J₂)   # effective trial-stress (von Mises stress)
    σᵗₑ = vonMises(σᵗ)
    σʸ = material.σ₀ + H * state.k # Previous yield limit
    # # println((ϵ[2, 2], state.ϵᵖ[2, 2], σᵗ[2, 2], σᵗₑ, σʸ))
    # println((state.ϵ[2, 2], Δϵ[2, 2], state.ϵᵖ[2, 2], vonMises(state.σ), σᵗₑ, σʸ))

    φᵗ  = σᵗₑ - σʸ # trial-value of the yield surface

    if φᵗ < 0.0 # elastic loading
        # println("elastic")
        # return σᵗ, material.Dᵉ, MaterialState(ϵ, state.ϵᵖ, σᵗ, state.k)
        return Δϵ, σᵗ, material.Dᵉ, MaterialState(state.ϵ + Δϵ, state.ϵᵖ, σᵗ, state.k)
    else # plastic loading
        h = H + 3G
        μ =  φᵗ / h   # plastic multiplier
        c1 = 1 - 3G * μ / σᵗₑ
        s = c1 * dev(σᵗ)           # updated deviatoric stress
        σ = s + vol(σᵗ)       # updated stress

        # compute algorithmic tangent stiffness ``D = \frac{\Delta \sigma }{\Delta \epsilon}``
        κ = H * (state.k + μ) # drag stress
        σₑ = material.σ₀ + κ  # updated yield surface
        Q(i,j,k,l) = Isymdev(i,j,k,l) - 3.0 / (2.0*σₑ^2) * s[i,j]*s[k,l]
        b = (3G*μ/σₑ) / (1.0 + 3G*μ/σₑ)
        Dtemp(i,j,k,l) = -2G*b * Q(i,j,k,l) - 9G^2 / (h*σₑ^2) * s[i,j]*s[k,l]
        D = material.Dᵉ + SymmetricTensor{4, 3}(Dtemp)

        # return new state
        Δϵᵖ = 3/2 * μ / σₑ * s # plastic strain
        ϵᵖ = state.ϵᵖ + Δϵᵖ    # plastic strain
        k = state.k + μ        # hardening variable
        # println(Δϵᵖ[2, 2])
        # println("PLASTIC")
        # return σ, D, MaterialState(ϵ, ϵᵖ, σ, k)
        return Δϵ, σ, D, MaterialState(state.ϵ + Δϵ, ϵᵖ, σ, k)
        # return Δϵ + Δϵᵖ, σ, D, MaterialState(state.ϵ + Δϵ + Δϵᵖ, ϵᵖ, σ, k) # has a shallower slope for some reason
    end
end

# function compute_stress_tangent(ϵ::SymmetricTensor{2, 3}, material::JCPlasticity, state::JohnsonCookState)
function compute_stress_tangent(Δϵ::SymmetricTensor{2, 3}, material::JCPlasticity, state::JohnsonCookState)
    # unpack some material parameters
    G = material.G

    # we use (•)ᵗ to denote *trial*-values
    # σᵗ = material.Dᵉ ⊡ (ϵ - state.ϵᵖ) # trial-stress
    # Δϵ = ϵ - state.ϵ

    # σᵗ = state.σ + material.Dᵉ ⊡ (ϵ - state.ϵ) # trial-stress
    # Δϵ = ϵ - state.ϵ
    σᵗ = state.σ + material.Dᵉ ⊡ Δϵ # trial-stress
    # σᵗ = material.Dᵉ ⊡ (state.ϵ + Δϵ) # trial-stress

    # sᵗ = dev(σᵗ)         # deviatoric part of trial-stress
    # J₂ = 0.5 * sᵗ ⊡ sᵗ   # second invariant of sᵗ
    # σᵗₑ = sqrt(3.0*J₂)   # effective trial-stress (von Mises stress)
    σᵗₑ = vonMises(σᵗ)
    σʸ = material.σ₀ # + H * state.k # Previous yield limit

    φᵗ  = σᵗₑ - σʸ # trial-value of the yield surface

    if φᵗ < 0.0 # elastic loading
        # return σᵗ, material.Dᵉ, JohnsonCookState(ϵ, state.ϵᵖ, σᵗ, state.k)
        return Δϵ, σᵗ, material.Dᵉ, JohnsonCookState(state.ϵ + Δϵ, state.ϵᵖ, σᵗ, state.k)
    else # plastic loading
        A   = 60.247247247247245
        B   = 520.9709709709709
        n   = 0.545
        C   = 0.022
        m   = 1.0
        Tr  = 298.0
        Tm  = 1356.0
        er0 = 1.0
        ϵ_dot = 10^-3
        θ = 295.0
        ϵ⁺  = ϵ_dot / er0
        θ⁺  = ( θ - Tr ) / ( Tm - Tr )
        h = 3G
        μ =  φᵗ / h   # plastic multiplier
        # # c1 = 1 - 3G * μ / σᵗₑ
        # # s = c1 * dev(σᵗ)           # updated deviatoric stress
        # # σ = s + vol(σᵗ)       # updated stress
        # ϵ = state.ϵ + Δϵ
        # # σ = SymmetricTensor{2, 3}(#= [ =#(  A .+ ( B .* (ϵ .^ n) )  ) .* (  ( 1. + C * log(ϵ⁺) ) * ( 1. - (θ⁺ ^ m) )  )#= ] =#)
        # σ = SymmetricTensor{2, 3}(#= [ =#(  A .+ ( B .* ( sign.(ϵ) .* ((abs.(ϵ) .^ n)) ) )  ) .* (  ( 1. + C * log(ϵ⁺) ) * ( 1. - ( sign(θ⁺) * ( abs(θ⁺) ^ m ) ) )  )#= ] =#)

        ϵₚᵗ = state.ϵᵖ + Δϵ
        σᵖ_eff = (#= [ =#(  A + ( B * ( vonMises(ϵₚᵗ) ^ n ) )  ) * (  ( 1.0 + C * log(ϵ⁺) ) * ( 1.0 - ( sign(θ⁺) * ( abs(θ⁺) ^ m ) ) )  )#= ] =#)
        # println(norm(ϵₚᵗ ./ norm(ϵₚᵗ)) == 1.)
        σ = SymmetricTensor{2, 3}( σᵖ_eff .* ( ϵₚᵗ ./ norm(ϵₚᵗ) ) )
        # n = dev(σᵗ) / vonMises(σᵗ)  # Normalized deviatoric trial stress direction
        # σ = SymmetricTensor{2, 3}( σᵖ_eff * (dev(σᵗ) / vonMises(σᵗ)) )  # Stress follows the plastic flow direction
        s = dev(σ)
        # # s = reduce(/, map(vonMises, [σ, state.σ])) * dev(state.σ)

        # compute algorithmic tangent stiffness ``D = \frac{\Delta \sigma }{\Delta \epsilon}``
        # κ = H * (state.k + μ) # drag stress
        σₑ = material.σ₀ # + κ  # updated yield surface
        Q(i,j,k,l) = Isymdev(i,j,k,l) - 3.0 / (2.0*σₑ^2) * s[i,j]*s[k,l]
        b = (3G*μ/σₑ) / (1.0 + 3G*μ/σₑ)
        Dtemp(i,j,k,l) = -2G*b * Q(i,j,k,l) - 9G^2 / (h*σₑ^2) * s[i,j]*s[k,l]
        D = material.Dᵉ + SymmetricTensor{4, 3}(Dtemp)
        # σₑ = material.σ₀ # updated yield surface
        # D = material.Dᵉ

        # return new state
        # Δϵᵖ = 3/2 * μ / σₑ * s # plastic strain
        Δϵᵖ = Δϵ # 3/2 * μ / σₑ * s # plastic strain
        ϵᵖ = state.ϵᵖ + Δϵᵖ    # plastic strain
        k = state.k + μ        # hardening variable
        # return σ, D, JohnsonCookState(ϵ, ϵᵖ, σ, k)
        return Δϵ, σ, D, JohnsonCookState(state.ϵ + Δϵ, ϵᵖ, σ, k)
        # return Δϵ, Δσ, D, JohnsonCookState(state.ϵ + Δϵ, ϵᵖ, Δσ, k)
        # # return Δϵ + Δϵᵖ, σ, D, JohnsonCookState(state.ϵ + Δϵ + Δϵᵖ, ϵᵖ, σ, k) # has a shallower slope for some reason
    end
end

function vector_projection(a::Vector{T}, b::Vector{T}) where T
    return (dot(a, b) / dot(b, b)) * b
end

do_verify = true
use_rve_or_E8 = true
if do_verify
    exp_data = CSV.read("data/verification/tension/Cu 500C_1_FullProcessed-true.csv", DataFrame; header=true, types=Float64)
    expdata_θ = 295.0
    expdata_ϵ̇ = 10^-3
    expdata_N = findlast(exp_data[!, "Strain"] .<= 0.35)
    expdata_strain = exp_data[!, "Strain"][1:expdata_N]
    expdata_displacement = (use_rve_or_E8 ? 1.0 : 32.0) .* expdata_strain
    expdata_time = expdata_strain ./ expdata_ϵ̇
    expdata_strainrate = (expdata_strain[begin+1:end] - expdata_strain[begin:end-1]) ./ (expdata_time[begin+1:end] - expdata_time[begin:end-1])
    pushfirst!(expdata_strainrate, 0.0)
    # expdata_displacement_temp = [first(expdata_displacement)]; for (i, x) in enumerate(expdata_displacement)
    #     if 1 < i < expdata_N
    #         append!(expdata_displacement_temp, [
    #             (expdata_displacement[i - 1] + x) / 2.0,
    #             x,
    #             (expdata_displacement[i + 1] + x) / 2.0,
    #         ])
    #     end
    # end; push!(expdata_displacement_temp, last(expdata_displacement)); expdata_displacement = expdata_displacement_temp
    expdata_stress = exp_data[!, "Stress (MPa)"][1:expdata_N]
    model_data = CSV.read("data/verification/tension/Cu 500C_1_JCdata.csv", DataFrame; header=true, types=Float64)
    modeldata_θ = 295.0
    modeldata_ϵ̇ = 10^-3
    modeldata_N = findlast(model_data[!, "strain-OFHC Cu"] .<= 0.35)
    modeldata_strain = model_data[!, "strain-OFHC Cu"][1:modeldata_N] # [vcat([1:5...], [10:10:modeldata_N...])]
    modeldata_displacement = (use_rve_or_E8 ? 1.0 : 32.0) .* modeldata_strain
    modeldata_time = modeldata_strain ./ modeldata_ϵ̇
    modeldata_strainrate = (modeldata_strain[begin+1:end] - modeldata_strain[begin:end-1]) ./ (modeldata_time[begin+1:end] - modeldata_time[begin:end-1])
    pushfirst!(modeldata_strainrate, 0.0)
    modeldata_stress = model_data[!, "VMstressOFHC Cu"][1:modeldata_N] # [vcat([1:5...], [10:10:modeldata_N...])]
end

function solve()
    simulationtime_0 = Dates.now()
    # define material parameters
    # E = 176.19474497681603e6    # [Pa]
    E = 95.28e3 * (!do_verify ? 1e6 : 1.0)
    # E = 120.0e3
    H = E/20                    # [Pa]
    ν = 0.343                     # [-]
    σ₀ = 62.996 * (!do_verify ? 1e6 : 1.0)               # [Pa]
    # material_Cu = J2Plasticity(E, ν, σ₀, H)
    material_Cu = JCPlasticity(E, ν, σ₀)
    # material_Fe = J2Plasticity(200e9, ν, 200e6, 200e9 / 20)
    material_Fe = JCPlasticity(200.0e3 * (!do_verify ? 1e6 : 1.0), ν, 200.0 * (!do_verify ? 1e6 : 1.0))
    material_endplate = material_Cu
    material_indentertip = material_Fe

    filename = do_verify ? (use_rve_or_E8 ? "rve-Cu" : "E8") : "plate-quarter_symmetry"
    # filename = "plate-quarter_symmetry-hexahedral"


    # define geometry
    if !do_verify
        ## general
        g                     ::Float64   = -9.80655          # [m/s^2]
        ## plate dimensions
        plate_width           ::Float64   = 0.0096            # [m]
        # plate_mass = 0.25                                           # [kg]
        ## ball dimenions
        # * measurements added 230110
        ball_diameter         ::Float64 = 0.0095              # [m]
        ball_radius           ::Float64 = ball_diameter/2     # [m]
        ball_rho              ::Float64 = 7667.3336           # [kg/m3] (https://www.salemball.com/440-stainless-steel-balls/)
        ball_nu               ::Float64 = 0.28
        ball_young            ::Float64 = 200e9               # [Pa]
        ball_yield            ::Float64 = 415e6               # [Pa]
        ball_mass             ::Float64 = 4pi/3*(             # [kg]
                ball_radius^3. * ball_rho)
        ball_ig               ::Float64 = 2ball_mass/5*(      # [kg-m2]
                ball_radius^2.)
        R = ball_radius
        mass = ball_mass/8
        # R = 6.35
        # mass = 26
        plate_diameter = 25.4/2/1000
        plate_thickness = 25.4/4/1000
        standoff_distance = 25.4/2/1000
    end


    # create geometry, dofs, and time-independent (constant) boundary conditions
    if do_verify
        if use_rve_or_E8
            L = 1.0
            grid = generate_grid(Hexahedron, (1, 1, 1), Vec{3}([0., 0., 0.]), Vec{3}([1., 1., 1.] .* L))
            addfacetset!(grid, "back_yz", x-> x[1] ≈ 0.)
            addfacetset!(grid, "back_xz", x-> x[2] ≈ 0.)
            addfacetset!(grid, "back_xy", x-> x[3] ≈ 0.)
            addfacetset!(grid, "front_xz", x-> x[2] ≈ L)
            # addcellset!(grid, "front_xz", x-> x[2] ≈ L)
            addnodeset!(grid, "rve", x->x[1] >= 0.)
            addnodeset!(grid, "front_xz", x->x[2] ≈ L)
            gridnodeset = getnodeset(grid, "rve")
            refshape = RefHexahedron
        else
            thickness = 3.          # Thickness of the specimen (B) in mm
            gauge_length = 32.      # Gauge length (L₀) in mm
            fillet_radius = 6.      # Fillet radius (R) in mm
            grip_width = 10.        # Grip width (nominal) in mm
            grip_length = 30.       # Grip length (in mm)
            gauge_width = 6.
            fillet_halfchord = √(fillet_radius^2 - (fillet_radius - (grip_width - gauge_width)/2)^2)
            # grid = togrid("data/mesh/E8.msh")
            grid = togrid("tensile_specimen.msh")
            L = maximum(x->get_node_coordinate(x)[1], getnodes(grid))
            # addcellset!(grid, "right", x->x[1] ≈ L)
            # addfacetset!(grid, "right", x-> x[1] ≈ L) # already exists
            addnodeset!(grid, "E8", 1:getnnodes(grid))
            addnodeset!(grid, "right", x->x[1] ≈ L)
            addnodeset!(grid, "gauge_left", x->x[1] ≈ grip_length + fillet_halfchord)
            addnodeset!(grid, "gauge_right", x->x[1] ≈ grip_length + fillet_halfchord + gauge_length)
            gridnodeset = getnodeset(grid, "E8")
            refshape = RefTetrahedron
        end
        gridncells = getncells(grid)
        gridnnodes = getnnodes(grid)
        interpolation = Lagrange{refshape, 1}()^3
    else
        grid = togrid(filename * ".msh")
        node_sphere_y = 7
        node_sphere_x = 8
        node_sphere_z = 9
        node_sphere_o = 10
        gridcellset_endplate = getcellset(grid, "end_plate")
        gridcellset_indentertip = getcellset(grid, "indenter_tip")
        # addfacetset!(grid, "cut_plane", x->x[3] ≈ plate_thickness + standoff_distance + R)
        addfacetset!(grid, "end_plate", x->(x[3] < plate_thickness + standoff_distance/2))
        addnodeset!(grid, "end_plate", x->(x[3] < plate_thickness + standoff_distance/2))
        addnodeset!(grid, "end_plate_cut_xz", x->((x[3] < plate_thickness + standoff_distance/2) && (x[2] ≈ 0.)))
        addnodeset!(grid, "end_plate_cut_yz", x->((x[3] < plate_thickness + standoff_distance/2) && (x[1] ≈ 0.)))
        addfacetset!(grid, "indenter_tip_all", x->(x[3] > plate_thickness + standoff_distance/2))
        addnodeset!(grid, "indenter_tip", x->(x[3] > plate_thickness + standoff_distance/2))
        gridnodeset_endplate = getnodeset(grid, "end_plate")
        gridnodeset_indentertip = getnodeset(grid, "indenter_tip")
        # (182, 306) [right, indenter_tip]
        # addnodeset!(grid, "left", x->x[3] ≈ 0)
        addnodeset!(grid, "right", x->x[3] ≈ plate_thickness)
        addnodeset!(grid, "indenter_tip_surface", x->abs(norm([0, 0, plate_thickness + standoff_distance + R] - x) - R) < R/24)
        addnodeset!(grid, "indenter_tip_cut_xy", x->x[3] ≈ plate_thickness + standoff_distance + R)
        gridnodeset_right = getnodeset(grid, "right")
        gridnodeset_indentertipsurface = getnodeset(grid, "indenter_tip_surface")
        gridncells = getncells(grid)
        gridncells_endplate = length(gridcellset_endplate)
        gridncells_indentertip = length(gridcellset_indentertip)
        gridnnodes = getnnodes(grid)
        gridnnodes_endplate = length(gridnodeset_endplate)
        gridnnodes_indentertip = length(gridnodeset_indentertip)
        refshape = RefTetrahedron
        interpolation = Lagrange{refshape, 1}()^3
    end

    ## degrees of freedom
    # dh = create_dofhandler(grid, interpolation) # JuaFEM helper function
    dh = DofHandler(grid)
    add!(dh, :u, interpolation) # add a displacement field with 3 components
    close!(dh)

    vectorstrafe_i(i) = vectorstrafe(grid, dh, i)

    ## constraints
    # dbcs = create_bc(dh, grid) # create Dirichlet boundary-conditions
    ch = ConstraintHandler(dh)
    if do_verify
        if use_rve_or_E8
            add!(ch, Dirichlet(:u, getfacetset(grid, "back_xz"), (x, t) -> [0.0], [2]))
            add!(ch, Dirichlet(:u, getfacetset(grid, "back_yz"), (x, t) -> [0.0], [1]))
            add!(ch, Dirichlet(:u, getfacetset(grid, "back_xy"), (x, t) -> [0.0], [3]))
        else
            add!(ch, Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3]))
        end
    else
        # // [20250124T1358] (JMA3): come back to this for proper constraint
        # [20250207T0954] (JMA3): pretty sure I fixed this last week (?) see "13. fixed displacement"
        # add!(ch, Dirichlet(:u, getfacetset(grid, "cylinder"), x -> [0., 0., 0.], [1, 2, 3]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "left"), x -> [0., 0., 0.], [1, 2, 3]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "end_plate_cut_xz"), x -> [0.], [2]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "end_plate_cut_yz"), x -> [0.], [1]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "indenter_tip_cut_xy"), x -> [0, 0], [1, 2]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "indenter_tip_cut_xz"), x -> [0.], [2]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "indenter_tip_cut_yz"), x -> [0.], [1]))
    end
    close!(ch)


    # define boundary conditions and time steps
    dN = 2
    n_timesteps = do_verify ? expdata_N : 1001
    dN != 1 ? (n_timesteps ÷= dN) : nothing
    displacement = do_verify ? collect(range(first(expdata_strain), last(expdata_strain); length=n_timesteps)) : nothing
    # V = -1
    V = do_verify ? expdata_ϵ̇ : -1.3477876909913353 # initial velocity onto Right plate in local reference frame
    time_0 = do_verify ? (first(expdata_strain) / V) : 0.0
    time_n = do_verify ? (last(expdata_strain) / V) : 0.1
    time_domain = collect(range(time_0, time_n; length=n_timesteps))
    # if time_0 != 0.
    #     n_timesteps += 1
    #     pushfirst!(time_domain,             0.0)
    #     pushfirst!(expdata_displacement,    0.0)
    #     pushfirst!(expdata_strain,          0.0)
    #     pushfirst!(expdata_stress,          0.0)
    # end
    if !do_verify
        time_domain = collect(time_domain)[begin:Int64((0.02 / (time_n / (n_timesteps - 1))) ÷ 1 + 1)]
        n_timesteps = length(time_domain)
    # else # if use_rve_or_E8
    #     # time_domain = collect(time_domain)[begin:Int64((0.025 / (time_n / (n_timesteps - 1))) ÷ 1 + 1)]
    #     # time_domain = collect(time_domain)[begin:Int64((2.0 / (time_n / (n_timesteps - 1))) ÷ 1 + 1)]
    #     time_domain = collect(time_domain)[begin:Int64((7.0 / (time_n / (n_timesteps - 1))) ÷ 1 + 1)]
    #     n_timesteps = length(time_domain)
    end
    VV = zeros(n_timesteps)
    ΔX = zeros(n_timesteps)


    # pre-allocate solution vectors, etc.
    # ! ndofs() ≡ dim * getnnodes(grid)
    n_dofs  = ndofs(dh)  # total number of dofs
    u       = zeros(n_dofs)  # solution vector
    Δu      = zeros(n_dofs)  # displacement correction
    u₀      = zeros(n_dofs)
    r       = zeros(n_dofs)   # residual
    K       = allocate_matrix(dh, ch); # tangent stiffness matrix


    # create material states. One array for each cell, where each element is an array of material-
    # states - one for each integration point
    ## quadratures
    # cellvalues_u, facetvalues_u = create_values(interpolation)
    # # cellvalues_u, cellvalues_eps, facetvalues = create_values(interpolation)
    # setup quadrature rules
    qr      = QuadratureRule{refshape}(2)
    facet_qr = FacetQuadratureRule{refshape}(refshape == RefHexahedron ? 2 : 3)
    # cell and facetvalues for u
    cellvalues = CellValues(qr, interpolation)
    facetvalues = FacetValues(facet_qr, interpolation)
    nu = getnbasefunctions(cellvalues)
    cell_basefuncs = getnbasefunctions(cellvalues)
    facet_basefuncs = getnbasefunctions(facetvalues)

    ## material states
    nqp_cell = getnquadpoints(cellvalues)
    nqp_facet = getnquadpoints(facetvalues)
    # TODO [20250214T0935] (JMA3) - figure out how to restart from previous calculation in case of interruptions
    # states_material = [MaterialState() for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_material_old = [MaterialState() for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_material = [JohnsonCookState() for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_material_old = [JohnsonCookState() for _ in 1:nqp_cell, _ in 1:gridncells]
    states_material = [[JohnsonCookState() for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:n_timesteps]
    # states_strain = [[SymmetricTensor{2, 3}(zeros((3, 3))) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells]
    states_strain = [[[SymmetricTensor{2, 3}(zeros((3, 3))) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:n_timesteps]
    # states_strain_old = [[SymmetricTensor{2, 3}(zeros((3, 3))) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_strain_rate = [[SymmetricTensor{2, 3}(zeros((3, 3))) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells]
    states_strain_rate = [[[SymmetricTensor{2, 3}(zeros((3, 3))) for _ in 1:cell_basefuncs] for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:n_timesteps]
    # states_energy_total = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_total_old = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_elastic = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_elastic_old = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_plastic = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    # states_energy_plastic_old = [0. for _ in 1:nqp_cell, _ in 1:gridncells]
    states_energy_total = [[0. for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:n_timesteps]
    states_energy_elastic = [[0. for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:n_timesteps]
    states_energy_plastic = [[0. for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:n_timesteps]
    # states_traction = [zeros((3, 3)) for _ in 1:nqp_cell, _ in 1:gridncells]
    states_traction = [[zeros((3, 3)) for _ in 1:nqp_cell, _ in 1:gridncells] for _ in 1:n_timesteps]
    # states_traction_old = [zeros((3, 3)) for _ in 1:nqp_cell, _ in 1:gridncells]


    # create export vectors results file
    u_max           = zeros(n_timesteps)
    strain_max      = zeros(n_timesteps)
    strainrate_max  = zeros(n_timesteps)
    stress_max      = zeros(n_timesteps)
    pvd             = paraview_collection(filename)


    gridnodes = deepcopy(getnodes(grid))
    gridnodes_view = @view gridnodes[[1:gridnnodes...]]
    if !do_verify
        gridnodes_endplate = gridnodes[[gridnodeset_endplate...]]
        gridnodes_endplate_view = @view gridnodes[[gridnodeset_endplate...]]
        gridnodes_indentertip = gridnodes[[gridnodeset_indentertip...]]
        gridnodes_indentertip_view = @view gridnodes[[gridnodeset_indentertip...]]
    end
    function update_rhs(grid, n)
        dx = get_node_coordinate(grid, n) - get_node_coordinate(gridnodes[node_sphere_o])
        Δx = (R * (dx / norm(dx))) - dx
        return Δx
    end
    function update_rhs(grid, n, v)
        px, py, pz = get_node_coordinate(gridnodes[n])
        cx, cy, cz = get_node_coordinate(gridnodes[node_sphere_o])
        vx, vy, vz = v
        # Calculate A, B, C for the quadratic equation
        A = vx^2 + vy^2 + vz^2
        B = 2 * ((px - cx) * vx + (py - cy) * vy + (pz - cz) * vz)
        C = (px - cx)^2 + (py - cy)^2 + (pz - cz)^2 - R^2

        # Discriminant of the quadratic equation
        discriminant = B^2 - 4A * C
        if discriminant < 0
            error("No real intersection: direction vector may be incorrect.")
        end

        # Solve for t using the quadratic formula
        t1 = (-B + sqrt(discriminant)) / (2A)
        t2 = (-B - sqrt(discriminant)) / (2A)

        # Choose the smallest positive t
        t_min = min(filter(x -> x > 0, [t1, t2])...)

        # Calculate the vector from P to the surface (Q)
        vector_to_surface = t_min * [vx, vy, vz]

        return vector_to_surface
    end
    function update_rhs_view(grid, n)
        return get_node_coordinate(gridnodes[node_sphere_z]) - get_node_coordinate(gridnodes[n])
    end


    # march through time
    NEWTON_TOL = 1 # 1 N
    NEWTON_M = 10
    println("\nStarting Netwon iterations:")
    # Δt = do_verify ? (time_domain[2] - time_domain[1]) : (time_n / (n_timesteps - 1))
    for (timestep, t) in enumerate(time_domain)
    # println((length(displacement), n_timesteps))
    # for (timestep, t) in zip(1 .+ (2 .* (1:((length(displacement) ÷ 2)) .- 1)), time_domain[begin:2:end])
        # @printf("\n Time step %d (t = %.6f s) @ %.4f m/s:\n", timestep - 1, t, V) # -1 to match ParaView
        VV[timestep] = V
        ΔV = 0.
        p = 0.
        # state_material = states_material
        # state_material_old = states_material_old
        # # println(first(state_material).σ == first(state_material_old).σ)
        state_material_old = (timestep > 1 ? states_material[timestep - 1] : states_material[1])
        state_material = (timestep > 1 ? states_material[timestep - 1] : states_material[1])
        # state_strain = states_strain
        # state_strain_old = states_strain_old
        state_strain_old = (timestep > 1 ? states_strain[timestep - 1] : states_strain[1])
        state_strain = (timestep > 1 ? states_strain[timestep - 1] : states_strain[1])
        state_strain_itr = deepcopy(state_strain)
        # state_energy_total = states_energy_total
        # state_energy_elastic = states_energy_elastic
        # state_energy_plastic = states_energy_plastic
        # state_energy_total_old = states_energy_total_old
        # state_energy_total_itr = deepcopy(state_energy_total_old)
        state_energy_total = states_energy_total[timestep]
        state_energy_elastic = states_energy_elastic[timestep]
        state_energy_plastic = states_energy_plastic[timestep]
        state_energy_total_old = (timestep > 1 ? states_energy_total[timestep - 1] : state_energy_total)
        state_energy_total_itr = deepcopy(state_energy_total_old)
        # state_traction = states_traction
        # state_traction_old = states_traction_old
        state_traction = states_traction[timestep]
        state_traction_old = (timestep > 1 ? states_traction[timestep - 1] : states_traction[1])
        # # (v)    only needed if ch[f(x, t)]     (v)
        # update!(ch, t) # evaluates the D-bndc at time t
        # # (^)                                   (^)
        apply!(u, ch)  # set the prescribed values in the solution vector
        if timestep > 1 # t > 0
            Δt = do_verify ? (t - time_domain[timestep - 1]) : (time_n / (n_timesteps - 1))
            ΔX[timestep] = do_verify ? displacement[timestep] : (ΔX[timestep - 1] + V * Δt)
            # ΔX[timestep] = ΔX[timestep - 1] + (do_verify ? (displacement[timestep] - displacement[timestep - 1]) : (V * Δt))
            if do_verify
                # @printf("\t\tExpected intermediate displacement increment of front_xz face: %.6f [m]\n", ΔX[timestep] - ΔX[timestep - 1])
            else
                @printf("\t\tExpected intermediate displacement increment of cut_xy face: %.6f [m]\n", V * Δt)
            end
            dbcs_t = ConstraintHandler(dh)
            if do_verify
                if use_rve_or_E8
                    add!(dbcs_t, Dirichlet(:u, getfacetset(grid, "front_xz"), (x, t) -> [ΔX[timestep]], [2]))
                    # add!(dbcs_t, Dirichlet(:u, getfacetset(grid, "front_xz"), (x, t) -> [-ν * ΔX[timestep], ΔX[timestep], -ν * ΔX[timestep]], [1, 2, 3]))
                else
                    add!(dbcs_t, Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> [ΔX[timestep]], [1]))
                end
            else
                add!(dbcs_t, Dirichlet(:u, [getnodeset(grid, "indenter_tip")...], (x, t) -> [ΔX[timestep]], [3]))
            end
            close!(dbcs_t)
            apply!(u, dbcs_t)  # set the prescribed values in the solution vector
        else
            Δt = do_verify ? (time_domain[2] - time_domain[1]) : (time_n / (n_timesteps - 1))
        end
        if do_verify
            if use_rve_or_E8
                @printf("\n Time step %d (t = %.6f s) @ %.6f mm:\n", timestep - 1, t, ΔX[timestep]) # -1 to match ParaView
                println_i_rve(i) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m]\n",
                    i, get_node_coordinate(gridnodes[i])[2], u[vectorstrafe_i(i)][2], ΔX[timestep])
                println_i_rve(i, trac, et, ee, ep, p, ΔV) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m] | T:%.3e [Pa] | Eᵗₘₐₓ:%.3e [J] | Eᵉₘₐₓ:%.3e [J] | Eᵖₘₐₓ:%.3e [J] | %%_diff(E):%.3f [%%] | ΔV:%.6f [m/s]\n",
                    i, get_node_coordinate(gridnodes[i])[2], u[vectorstrafe_i(i)][2], ΔX[timestep], maximum(abs, map(norm, trac)), maximum(abs, et), maximum(abs, ee), maximum(abs, ep), p, ΔV)
            else
                @printf("\n Time step %d (t = %.6f s) @ %.6f mm/mm, %.6f mm:\n", timestep - 1, t, expdata_strain[timestep], ΔX[timestep]) # -1 to match ParaView
                println_i_E8(i) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m]\n",
                    i, get_node_coordinate(gridnodes[i])[1], u[vectorstrafe_i(i)][1], ΔX[timestep])
                println_i_E8(i, trac, et, ee, ep, p, ΔV) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m] | T:%.3e [Pa] | Eᵗₘₐₓ:%.3e [J] | Eᵉₘₐₓ:%.3e [J] | Eᵖₘₐₓ:%.3e [J] | %%_diff(E):%.3f [%%] | ΔV:%.6f [m/s]\n",
                    i, get_node_coordinate(gridnodes[i])[1], u[vectorstrafe_i(i)][1], ΔX[timestep], maximum(abs, map(norm, trac)), maximum(abs, et), maximum(abs, ee), maximum(abs, ep), p, ΔV)
            end
        else
            @printf("\n Time step %d (t = %.6f s) @ %.4f m/s:\n", timestep - 1, t, V) # -1 to match ParaView
            println_i(i, uprhs) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m] | Δx:%.6f [m]\n",
                i, get_node_coordinate(gridnodes[i])[3], u[vectorstrafe_i(i)][3], ΔX[timestep], uprhs(grid, i)[3])
            println_i(i, trac, et, ee, ep, p, ΔV) = @printf("\ti:%05d | n₃:%.6f [m] | u₃:%.6f [m] | ΔX:%.6f [m] | Δx:%.6f [m] | T:%.3e [Pa] | Eᵗₘₐₓ:%.3e [J] | Eᵉₘₐₓ:%.3e [J] | Eᵖₘₐₓ:%.3e [J] | %%_diff(E):%.3f [%%] | ΔV:%.6f [m/s]\n",
                i, get_node_coordinate(gridnodes[i])[3], u[vectorstrafe_i(i)][3], ΔX[timestep], update_rhs(grid, i)[3], maximum(abs, map(norm, trac)), maximum(abs, et), maximum(abs, ee), maximum(abs, ep), p, ΔV)
        end
        # (-sqrt(state_plasticenergy[q_point] / mass) * V) - V
        if do_verify
            # if use_rve_or_E8
            #     i = 8; println_i_rve(i)
            # else
            #     i = 23; println_i_E8(i)
            # end
        else
            if 80 < timestep < 140
                i = node_sphere_o; println_i(i, update_rhs)
                i = 5; println_i(i, update_rhs)
            end
        end


        if do_verify
            # if use_rve_or_E8
            #     for (i, j, n) in zip(getnodeset(grid, "front_xz"), eachindex(getnodeset(grid, "front_xz")), gridnodes[[getnodeset(grid, "front_xz")...]])
            #         gridnodes_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)]))
            #     end
            #     @printf("\t\tIntermediate displacement of front_xz face: %.4f [m]\n", average([u[vectorstrafe_i(n)][3] for n in getnodeset(grid, "front_xz")]))
            # else
            #     for (i, j, n) in zip(getnodeset(grid, "right"), eachindex(getnodeset(grid, "right")), gridnodes[[getnodeset(grid, "right")...]])
            #         gridnodes_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)]))
            #     end
            #     @printf("\t\tIntermediate displacement of right face: %.4f [m]\n", average([u[vectorstrafe_i(n)][1] for n in getnodeset(grid, "right")]))
            # end
        else
            # * [20250121T1533] (JMA3) - keep the `for`-loop below for properly strides through ⃗u to update node coordinates
            for (i, j, n) in zip(gridnodeset_endplate, 1:gridnnodes_endplate, gridnodes_endplate)
                gridnodes_endplate_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)]))
            end
            for (i, j, n) in zip(gridnodeset_indentertip, 1:gridnnodes_indentertip, gridnodes_indentertip)
                gridnodes_indentertip_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)]))
            end
            @printf("\t\tIntermediate displacement of cut_xy face: %.4f [m]\n", average([u[vectorstrafe_i(n)][3] for n in getnodeset(grid, "indenter_tip_cut_xy")]))
            itxnodes_endplate = Int64[]
            for (j, i, n) in zip(eachindex(gridnodeset_endplate), gridnodeset_endplate, #= TODO maybe this should be gridnodes_endplate_view? =# gridnodes_endplate_view)
                if norm(get_node_coordinate(n) - get_node_coordinate(gridnodes[node_sphere_o])) < R
                    push!(itxnodes_endplate, i)
                end
            end
            unique!(itxnodes_endplate)
            # intersect!(itxnodes_surface_endplate, getnodeset(grid, "right"))
            if !isempty(itxnodes_endplate)
                intersecting_coords = length(itxnodes_endplate)
                println("\t\titx coords: $intersecting_coords")
                addnodeset!(grid, "itxnodes-endplate-$timestep", itxnodes_endplate)
                itxnodes_facets_endplate = FacetIndex[]
                itxnodes_cells_endplate = Int64[]
                for fc in getfacetset(grid, "right")
                    if !isempty(intersect(getcells(grid)[fc[1]].nodes, itxnodes_endplate))
                        push!(itxnodes_facets_endplate, fc)
                        push!(itxnodes_cells_endplate, fc[1])
                    end
                end
                addfacetset!(grid, "itxnodes-endplate-$timestep", itxnodes_facets_endplate)
                addcellset!(grid, "itxnodes-endplate-$timestep", itxnodes_cells_endplate)
                # for n in getnodeset(grid, "itxnodes-endplate-$timestep")
                #     u[vectorstrafe_i(n)] += update_rhs(grid, n) # - Δu[vectorstrafe_i(n)]
                # end
                itxnodes_surface_endplate = union(Int64[], map(s->intersect(itxnodes_endplate, getnodeset(grid, s)), ["right", "end_plate_cut_xz", "end_plate_cut_yz"])...)
                intersecting_coords_surface = length(itxnodes_surface_endplate)
                println("\t\titx coords (surface): $intersecting_coords_surface")
                addnodeset!(grid, "itxnodes-endplate_right-$timestep", itxnodes_surface_endplate)
                itxnodes_surface_facets_endplate = FacetIndex[]
                itxnodes_surface_cells_endplate = Int64[]
                for fc in getfacetset(grid, "right")
                    if !isempty(intersect(getcells(grid)[fc[1]].nodes, itxnodes_surface_endplate))
                        push!(itxnodes_surface_facets_endplate, fc)
                        push!(itxnodes_surface_cells_endplate, fc[1])
                    end
                end
                addfacetset!(grid, "itxnodes-endplate_right-$timestep", itxnodes_surface_facets_endplate)
                addcellset!(grid, "itxnodes-endplate_right-$timestep", itxnodes_surface_cells_endplate)
                for n in getnodeset(grid, "itxnodes-endplate_right-$timestep")
                    u[vectorstrafe_i(n)] += update_rhs(grid, n) # - Δu[vectorstrafe_i(n)]
                end

                # # println((get_node_coordinate(gridnodes[node_sphere_o]), get_node_coordinate(gridnodes[getnodeset(grid, "itxnodes-endplate-$timestep.$newton_itr")[end]])))
                # for i in getnodeset(grid, "itxnodes-endplate-$timestep")
                #     println((get_node_coordinate(gridnodes[node_sphere_o]), get_node_coordinate(gridnodes[i])))
                # end

                # # return_itrvalue = any([(norm(get_node_coordinate(grid, n) - get_node_coordinate(gridnodes[node_sphere_o])) < R) for n in getnodeset(grid, "itxnodes-endplate_right-$timestep")])
                # # return_itr = 0
                # # while return_itrvalue
                # #     # * [20250121T1533] (JMA3) - keep the `for`-loop below for properly strides through ⃗u to update node coordinates
                # #     # for (i, j, n) in zip(gridnodeset_endplate, 1:gridnnodes_endplate, gridnodes_endplate)
                # #     #     gridnodes_endplate_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)]))
                # #     # end
                # #     # Δx = 0.
                # #     # for (i, j, n) in zip(gridnodeset_indentertip, 1:gridnnodes_indentertip, gridnodes_indentertip)
                # #     #     uprhs = update_rhs_view(grid, i)[3]
                # #     #     Δx = (Δx, uprhs)[findmax(abs, [Δx, uprhs])[2]]
                # #     # end
                # #     Δx = get_node_coordinate(gridnodes[node_sphere_z])[3] - get_node_coordinate(gridnodes[5])[3]
                # #     for (i, j, n) in zip(gridnodeset_indentertip, 1:gridnnodes_indentertip, gridnodes_indentertip)
                # #         # # gridnodes_indentertip_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)]))
                # #         gridnodes_indentertip_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)] + [0., 0., Δx]))
                # #     end
                # #     ΔX[timestep] = u[vectorstrafe_i(node_sphere_o)][3] - Δx
                # #     for i in gridnodeset_indentertip
                # #         u[vectorstrafe_i(i)] -= [0., 0., Δx]
                # #     end
                # #     return_itrvalue = any([(norm(get_node_coordinate(grid, n) - get_node_coordinate(gridnodes[node_sphere_o])) < R) for n in getnodeset(grid, "itxnodes-endplate_right-$timestep")])
                # #     return_itr += 1
                # # end
                # # # apply!(u, ch)  # set the prescribed values in the solution vector
                # # println("\tNumber of iterations to put ball back: ", return_itr)
                # Δx = get_node_coordinate(gridnodes[node_sphere_z])[3] - get_node_coordinate(gridnodes[5])[3]
                # # Δx = maximum(n->abs(update_rhs(grid, n, [0., 0., 1.])[3]), [getnodeset(grid, "itxnodes-endplate-$timestep")...])
                # for (i, j, n) in zip(gridnodeset_indentertip, 1:gridnnodes_indentertip, gridnodes_indentertip)
                #     gridnodes_indentertip_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)] + [0., 0., Δx]))
                #     # gridnodes_indentertip_view[j] = Node(Ferrite.Vec{3}(get_node_coordinate(n) + u[vectorstrafe_i(i)] + [0., 0., Δx]))
                # end
                # ΔX[timestep] = u[vectorstrafe_i(node_sphere_o)][3] - Δx
                # # ΔX[timestep] = u[vectorstrafe_i(node_sphere_o)][3] + Δx
                # for i in gridnodeset_indentertip
                #     u[vectorstrafe_i(i)] -= [0., 0., Δx]
                #     # u[vectorstrafe_i(i)] += [0., 0., Δx]
                # end
                # if do_verify
                #     if use_rve_or_E8
                #         i = 8; println_i_rve(i)
                #     else
                #         i = 23; println_i_E8(i)
                #     end
                # else
                #     if 80 < timestep < 140
                #         i = node_sphere_o; println_i(i, update_rhs)
                #     end
                # end
            end
        end


        # newton-raphson loop
        newton_itr = 1
        while newton_itr <= NEWTON_M
            if do_verify
                # Tangent and residual contribution from the cells (volume integral)
                # doassemble!(K, r, cellvalues, dh, material, u, states, states_old);
                assembler = start_assemble(K, r)
                # assembler = start_assemble(K, r; fillzero=false)
                if t > 0.
                    re = zeros(nu)     # element residual vector
                    ke = zeros(nu, nu) # element tangent matrix
                    for cell in CellIterator(dh)
                        fill!(ke, 0)
                        fill!(re, 0)
                        eldofs = celldofs(cell)
                        ue = u[eldofs]
                        # re = @view r[eldofs]
                        # ke = @view K[eldofs, eldofs]
                        state = @view state_material[:, cell.cellid]
                        state_old = @view state_material_old[:, cell.cellid]
                        strain = @view state_strain[:, cell.cellid]
                        strain_itr = @view state_strain_itr[:, cell.cellid]
                        strain_old = @view state_strain_old[:, cell.cellid]
                        material = material_Cu
                        stateenergy_total = @view state_energy_total[:, cell.cellid]
                        stateenergy_elastic = @view state_energy_elastic[:, cell.cellid]
                        stateenergy_plastic = @view state_energy_plastic[:, cell.cellid]
                        statetraction = @view state_traction[:, cell.cellid]
                        fill!(stateenergy_total,    0.0)
                        fill!(stateenergy_elastic,  0.0)
                        fill!(stateenergy_plastic,  0.0)
                        # fill!(statetraction,        0.0)
                        Ferrite.reinit!(cellvalues, cell)
                        for q_point in 1:nqp_cell
                            # for each integration point, compute stress and material stiffness
                            # ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
                            # σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state_old[q_point])
                            Δϵ = function_symmetric_gradient(cellvalues, q_point, ue) - state[q_point].ϵ # Total strain
                            σ_prev, ϵ_prev, ϵᵖ_prev = state[q_point].σ, state[q_point].ϵ, state[q_point].ϵᵖ
                            # if q_point == 1
                            #     # println(size(function_symmetric_gradient(cellvalues, q_point, ue)))
                            #     println(size(shape_gradient(cellvalues, q_point, 1)))
                            # end
                            Δϵ, σ, D, state[q_point] = compute_stress_tangent(Δϵ, material, state[q_point])
                            ϵ = state[q_point].ϵ
                            # σ_prev, ϵ_prev, ϵᵖ_prev = state[q_point].σ, state[q_point].ϵ, state[q_point].ϵᵖ
                            # ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
                            # σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state[q_point])
                            # Δϵ = ϵ - ϵ_prev
                            dΩ = getdetJdV(cellvalues, q_point)
                            # stateenergy_total[q_point]      = (transpose(ϵ) ⊡ σ) * dΩ
                            # stateenergy_elastic[q_point]    = (transpose(ϵ - state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            # stateenergy_plastic[q_point]    = (transpose(state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            # statetraction[q_point]          = σ # - σ_prev
                            stateenergy_total[q_point]     += (transpose(Δϵ) ⊡ (σ - σ_prev)) * dΩ
                            stateenergy_elastic[q_point]   += (transpose((ϵ - state[q_point].ϵᵖ) - (ϵ_prev - ϵᵖ_prev)) ⊡ (σ - σ_prev)) * dΩ
                            stateenergy_plastic[q_point]   += (transpose(state[q_point].ϵᵖ - ϵᵖ_prev) ⊡ (σ - σ_prev)) * dΩ
                            statetraction[q_point]         += σ # - σ_prev
                            for i in 1:cell_basefuncs
                                δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
                                # re[i] += (δϵ ⊡ σ) * dΩ # TODO this needs to be uncommented! # add internal force to residual
                                re[i] += ((δϵ - strain[q_point][i]) ⊡ (σ - σ_prev)) * dΩ
                                for j in 1:i # loop only over lower half
                                    Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                                    # ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
                                    ke[i, j] += (δϵ - strain[q_point][i]) ⊡ D ⊡ (Δϵ - strain[q_point][j]) * dΩ
                                end
                                strain[q_point][i] = δϵ
                            end
                        end
                        symmetrize_lower!(ke)
                        assemble!(assembler, eldofs, ke, re)
                    end

                    # Residual contribution from the Neumann boundary (surface integral)
                    # doassemble_neumann!(r, dh, getfacetset(grid, "right"), facetvalues, traction)
                    # n_basefuncs = getnbasefunctions(facetvalues)
                    rf = zeros(facet_basefuncs)                      # element residual vector
                    for fc in FacetIterator(dh, getfacetset(grid, use_rve_or_E8 ? "front_xz" : "right"))
                        # Add traction as a negative contribution to the element residual `re`:
                        reinit!(facetvalues, fc)
                        fill!(rf, 0)
                        # rf = @view r[celldofs(fc)]
                        # rf .= 0.0
                        # uf = @view u[celldofs(fc)]
                        stateenergy_total = state_energy_total[:, fc.cc.cellid]
                        stateenergy_total_itr = state_energy_total_itr[:, fc.cc.cellid]
                        stateenergy_total_old = state_energy_total_old[:, fc.cc.cellid]
                        statetraction = state_traction[:, fc.cc.cellid]
                        for q_point in 1:nqp_facet
                            traction_u = (stateenergy_total[q_point] - stateenergy_total_itr[q_point]) * getnormal(facetvalues, q_point)
                            # traction_u = stateenergy_total[q_point] * getnormal(facetvalues, q_point)
                            # traction_u = stateenergy_total[q_point] * Ferrite.Vec{3}([1.0, 0.0, 0.0])
                            # traction_u = (stateenergy_total[q_point] - stateenergy_total_itr[q_point]) * Ferrite.Vec{3}([0.0, 1.0, 0.0])
                            # traction_u = (stateenergy_total[q_point] - stateenergy_total_old[q_point]) * Ferrite.Vec{3}([1.0, 0.0, 0.0])
                            # traction_u = statetraction[q_point] * getnormal(facetvalues, q_point)
                            # traction_u = statetraction[q_point] * Ferrite.Vec{3}([0.0, 1.0, 0.0])
                            # traction_u = Ferrite.Vec{3}(3000e3 .* [1.0, 0.0, 0.0])
                            dΓ = getdetJdV(facetvalues, q_point)
                            # rf[q_point] -= (uf[3(q_point - 1) .+ [1, 2, 3]] ⋅ traction_u) * dΓ
                            for i in 1:facet_basefuncs
                                δu = shape_value(facetvalues, q_point, i)
                                rf[i] -= (δu ⋅ traction_u) * dΓ
                            end
                        end
                        assemble!(r, celldofs(fc), rf)
                    end


                    state_strain_itr       .= state_strain
                    state_energy_total_itr .= state_energy_total
                end
                # if use_rve_or_E8
                #     i = 8; println_i_rve(i, state_traction, state_energy_total, state_energy_elastic, state_energy_plastic, p, ΔV)
                # else
                #     i = 23; println_i_E8(i, state_traction, state_energy_total, state_energy_elastic, state_energy_plastic, p, ΔV)
                # end
            else
                if !isempty(itxnodes_endplate)
                    # tangent and residual contribution from the cells (volume integral)
                    assembler = start_assemble(K, r)
                    traction = [Ferrite.Vec{3}((0., 0., 0.)) for _ in 1:nqp_cell, _ in 1:gridncells]

                    re = zeros(nu)     # element residual vector
                    ke = zeros(nu, nu) # element tangent matrix
                    for cell in CellIterator(dh, gridcellset_endplate)
                        fill!(ke, 0)
                        fill!(re, 0)
                        eldofs = celldofs(cell)
                        ue = u[eldofs]
                        # ue = u[eldofs] - u₀[eldofs]
                        state = @view state_material[:, cell.cellid]
                        state_old = @view state_material_old[:, cell.cellid]
                        strain = @view state_strain[:, cell.cellid]
                        material = material_endplate
                        stateenergy_total = @view state_energy_total[:, cell.cellid]
                        stateenergy_elastic = @view state_energy_elastic[:, cell.cellid]
                        stateenergy_plastic = @view state_energy_plastic[:, cell.cellid]
                        fill!(stateenergy_total,    0.0)
                        fill!(stateenergy_elastic,  0.0)
                        fill!(stateenergy_plastic,  0.0)
                        statetraction = @view state_traction[:, cell.cellid]
                        Ferrite.reinit!(cellvalues, cell)
                        for q_point in 1:nqp_cell
                            Ferrite.reinit!(facetvalues, cell, q_point)
                            # For each integration point, compute stress and material stiffness
                            # ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # , dof_range_endplate) # Total strain
                            # σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state_old[q_point])
                            # dΩ = getdetJdV(cellvalues, q_point)
                            # stateenergy_total[q_point] = 0.5 * (transpose(ϵ) ⊡ σ) * dΩ
                            # stateenergy_elastic[q_point] = 0.5 * (transpose(ϵ - state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            # stateenergy_plastic[q_point] = 0.5 * (transpose(state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            # # stateenergy_total[q_point] = 0.5 * (transpose(ϵ) ⊡ D ⊡ ϵ) * dΩ
                            # # stateenergy_elastic[q_point] = 0.5 * (transpose(ϵ - state[q_point].ϵᵖ) ⊡ D ⊡ (ϵ - state[q_point].ϵᵖ)) * dΩ
                            # # stateenergy_plastic[q_point] = 0.5 * (transpose(state[q_point].ϵᵖ) ⊡ D ⊡ state[q_point].ϵᵖ) * dΩ
                            # # traction[q_point, cell.cellid] += (0.5 * transpose(ϵ - state[q_point].ϵᵖ) ⊡ material.Dᵉ ⊡ (ϵ - state[q_point].ϵᵖ)) * getnormal(facetvalues, q_point) # * dΩ
                            # traction[q_point, cell.cellid] += ((transpose(ϵ - state[q_point].ϵᵖ) ⊡ σ)) * getnormal(facetvalues, q_point) # * dΩ
                            # # traction[q_point, cell.cellid] += Ferrite.Vec{3}(vector_projection([((0.5 * (transpose(ϵ - state[q_point].ϵᵖ) ⊡ σ)) * getnormal(facetvalues, q_point))...], sign(V) .* [0., 0., 1.]))
                            σ_prev, ϵ_prev, ϵᵖ_prev = state[q_point].σ, state[q_point].ϵ, state[q_point].ϵᵖ
                            ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
                            σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state[q_point])
                            Δϵ = ϵ - ϵ_prev
                            dΩ = getdetJdV(cellvalues, q_point)
                            stateenergy_total[q_point]      = (transpose(ϵ) ⊡ σ) * dΩ
                            stateenergy_elastic[q_point]    = (transpose(ϵ - state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            stateenergy_plastic[q_point]    = (transpose(state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            statetraction[q_point]          = σ # - σ_prev
                            # stateenergy_total[q_point]     += (transpose(Δϵ) ⊡ (σ - σ_prev)) * dΩ
                            # stateenergy_elastic[q_point]   += (transpose((state[q_point].ϵ - state[q_point].ϵᵖ) - (ϵ_prev - ϵᵖ_prev)) ⊡ (σ - σ_prev)) * dΩ
                            # stateenergy_plastic[q_point]   += (transpose(state[q_point].ϵᵖ - ϵᵖ_prev) ⊡ (σ - σ_prev)) * dΩ
                            # statetraction[q_point]         += σ - σ_prev
                            for i in 1:cell_basefuncs
                                δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
                                re[i] += (δϵ ⊡ σ) * dΩ # add internal force to residual
                                # re[i] += ((δϵ - strain[q_point][i]) ⊡ (σ - σ_prev)) * dΩ
                                for j in 1:i # loop only over lower half
                                    Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                                    ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
                                    # ke[i, j] += (δϵ - strain[q_point][i]) ⊡ D ⊡ (Δϵ - strain[q_point][j]) * dΩ
                                end
                                strain[q_point][i] = δϵ
                            end
                        end
                        symmetrize_lower!(ke)
                        assemble!(assembler, eldofs, ke, re)
                    end
                    re = zeros(nu)     # element residual vector
                    ke = zeros(nu, nu) # element tangent matrix
                    for cell in CellIterator(dh, gridcellset_indentertip)
                        fill!(ke, 0)
                        fill!(re, 0)
                        eldofs = celldofs(cell)
                        ue = u[eldofs] # [dof_range_indentertip]] # - u₀[eldofs[dof_range_indentertip]]
                        state = @view state_material[:, cell.cellid]
                        state_old = @view state_material_old[:, cell.cellid]
                        strain = @view state_strain[:, cell.cellid]
                        material = material_indentertip
                        Ferrite.reinit!(cellvalues, cell)
                        for q_point in 1:nqp_cell
                            Ferrite.reinit!(facetvalues, cell, q_point)
                            # For each integration point, compute stress and material stiffness
                            # ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # , dof_range_indentertip) # Total strain
                            # σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state_old[q_point])
                            # dΩ = getdetJdV(cellvalues, q_point)
                            # # traction[q_point, cell.cellid] += σ ⋅ getnormal(facetvalues, q_point) # + Ferrite.Vec{3}([0., 0., V + ΔV] .* (mass / Δt))
                            # # traction[q_point, cell.cellid] += Ferrite.Vec{3}(vector_projection([(σ ⋅ getnormal(facetvalues, q_point))...], V .* [0, 0, 1])) # + Ferrite.Vec{3}([0., 0., V + ΔV] .* (mass / Δt))
                            # # traction[q_point, cell.cellid] += Ferrite.Vec{3}([0., 0., V] .* (mass / Δt)) / dΩ
                            # traction[q_point, cell.cellid] += Ferrite.Vec{3}(vector_projection((sign(V) * mass) .* ((V .* [0, 0, 1]) .^ 2.), [getnormal(facetvalues, q_point)...]))
                            # # traction[q_point, cell.cellid] += Ferrite.Vec{3}(((sign(V) * mass) .* ((V .* [0, 0, 1]) .^ 2.))) ⋅ getnormal(facetvalues, q_point)
                            # # traction[q_point, cell.cellid] += 0.5 * Ferrite.Vec{3}(mass .* ((V .* [0, 0, 1]) .^ 2.)) * sign(V)
                            σ_prev, ϵ_prev, ϵᵖ_prev = state[q_point].σ, state[q_point].ϵ, state[q_point].ϵᵖ
                            ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
                            σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state[q_point])
                            Δϵ = ϵ - ϵ_prev
                            dΩ = getdetJdV(cellvalues, q_point)
                            # # stateenergy_total[q_point]      = (transpose(ϵ) ⊡ σ) * dΩ
                            # # stateenergy_elastic[q_point]    = (transpose(ϵ - state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            # # stateenergy_plastic[q_point]    = (transpose(state[q_point].ϵᵖ) ⊡ σ) * dΩ
                            # # statetraction[q_point]          = σ # - σ_prev
                            # stateenergy_total[q_point]     += (transpose(Δϵ) ⊡ (σ - σ_prev)) * dΩ
                            # stateenergy_elastic[q_point]   += (transpose((state[q_point].ϵ - state[q_point].ϵᵖ) - (ϵ_prev - ϵᵖ_prev)) ⊡ (σ - σ_prev)) * dΩ
                            # stateenergy_plastic[q_point]   += (transpose(state[q_point].ϵᵖ - ϵᵖ_prev) ⊡ (σ - σ_prev)) * dΩ
                            # statetraction[q_point]         += σ - σ_prev
                            for i in 1:cell_basefuncs
                                δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
                                # traction[q_point, cell.cellid] += 0.5 * ((δϵ ⊡ σ) * getnormal(facetvalues, q_point) - Ferrite.Vec{3}(vector_projection(mass .* ((V .* [0, 0, 1]) .^ 2.), [getnormal(facetvalues, q_point)...])))
                                re[i] += (-δϵ ⊡ σ) * dΩ # add internal force to residual
                                # re[i] += (-(δϵ - strain[q_point][i]) ⊡ (σ - σ_prev)) * dΩ
                                # traction[q_point, cell.cellid] -= (δϵ ⊡ σ) * getnormal(facetvalues, q_point)
                                for j in 1:i # loop only over lower half
                                    Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                                    ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
                                    # ke[i, j] += (δϵ - strain[q_point][i]) ⊡ D ⊡ (Δϵ - strain[q_point][j]) * dΩ
                                end
                                strain[q_point][i] = δϵ
                            end
                        end
                        symmetrize_lower!(ke)
                        assemble!(assembler, eldofs, ke, re)
                    end

                    # # residual contribution from the Neumann boundary (surface integral)
                    # # doassemble_neumann!(r, dh, getfacetset(grid, "left"), facetvalues_u, traction)
                    # rf = zeros(facet_basefuncs)                      # element residual vector
                    # for fc in FacetIterator(dh, getfacetset(grid, "itxnodes-endplate_right-$timestep"))
                    #     # add traction as a negative contribution to the element residual `re`:
                    #     Ferrite.reinit!(facetvalues, fc)
                    #     fill!(rf, 0)
                    #     for q_point in 1:nqp_facet
                    #         dΓ = getdetJdV(facetvalues, q_point)
                    #         traction_u = vector_projection([traction[q_point, fc.cc.cellid]...], sign(V) .* [0., 0., 1.])
                    #         # traction_u = traction[q_point, fc.cc.cellid]
                    #         for i in 1:facet_basefuncs
                    #             δu = shape_value(facetvalues, q_point, i)
                    #             rf[i] -= (δu ⋅ traction_u) * dΓ
                    #         end
                    #     end
                    #     assemble!(r, celldofs(fc), rf)
                    # end


                    # V = 1
                    ΔV₁ = 0.
                    V₂ = 0.
                    for cell in CellIterator(grid, getcellset(grid, "end_plate"))
                        Ferrite.reinit!(cellvalues, cell)
                        cell_id = cellid(cell)
                        # eldofs = celldofs(cell)
                        stateenergy_total = @view state_energy_total[:, cell_id]
                        # stateenergy_elastic = @view state_energy_elastic[:, cell_id]
                        # stateenergy_plastic = @view state_energy_plastic[:, cell_id]
                        ΔV₂ = 0.
                        for q_point in 1:nqp_cell
                            ΔV₂ += -sign(V) * sign(stateenergy_total[q_point]) * sqrt(abs(stateenergy_total[q_point]) / mass)
                            # ΔV₂ += sign(V) * sqrt(stateenergy_plastic[q_point] / mass)
                        end
                        # V₂ += V + ΔV₂ / 4
                        ΔV₁ += ΔV₂ / 4
                    end
                    # V₂ /= length(getcellset(grid, "end_plate"))
                    # V = V₂
                    ΔV₁ /= length(getcellset(grid, "end_plate"))
                    V += ΔV₁
                    ΔV₁ = 0.
                    # if abs(V) < 0.001 || percentdifference(map(x->findmax(abs, x)[1], [states_elasticenergy, states_elasticenergy_old])...) < 0.01
                    # if percentdifference(map(x->findmax(abs, x)[1], [state_energy_total, state_energy_total_old])...) < 1.
                    # p = if timestep > 4
                    #     maxenergystates = map(x->findmax(abs, x)[1], states_energy_total[timestep-4:timestep])
                    #     percdiff = zeros(length(maxenergystates) - 1)
                    #     for (i, mes1, mes2) in zip(eachindex(maxenergystates) - 1, maxenergystates[begin:end-1], maxenergystates[begin+1:end])
                    #         percdiff[i] = percentdifference(mes1, mes2)
                    #     end
                    #     norm(percdiff)
                    # elseif timestep > 3
                    #     maxenergystates = map(x->findmax(abs, x)[1], states_energy_total[timestep-3:timestep])
                    #     percdiff = zeros(length(maxenergystates) - 1)
                    #     for (i, mes1, mes2) in zip(eachindex(maxenergystates) - 1, maxenergystates[begin:end-1], maxenergystates[begin+1:end])
                    #         percdiff[i] = percentdifference(mes1, mes2)
                    #     end
                    #     norm(percdiff)
                    # elseif timestep > 2
                    p = if timestep > 2
                        maxenergystates = map(x->maximum(abs, x), states_energy_total[timestep-2:timestep])
                        percdiff = zeros(length(maxenergystates) - 1)
                        for (i, mes1, mes2) in zip(eachindex(maxenergystates), maxenergystates[begin:end-1], maxenergystates[begin+1:end])
                            percdiff[i] = percentdifference(mes1, mes2)
                        end
                        norm(percdiff)
                    elseif timestep > 1
                        percentdifference(map(x->maximum(abs, x), [state_energy_total, state_energy_total_old])...)
                    else
                        error("Could not determine elastic energy state at t[", timestep, "] = ", t, " s")
                    end
                    if p < 5
                        for cell in CellIterator(grid, getcellset(grid, "end_plate"))
                            Ferrite.reinit!(cellvalues, cell)
                            cell_id = cellid(cell)
                            # eldofs = celldofs(cell)
                            stateenergy_total = @view state_energy_total[:, cell_id]
                            stateenergy_elastic = @view state_energy_elastic[:, cell_id]
                            stateenergy_plastic = @view state_energy_plastic[:, cell_id]
                            ΔV₂ = 0.
                            for q_point in 1:nqp_cell
                                Ferrite.reinit!(facetvalues, cell, q_point)
                                # ΔV₂ += -sign(V) * sign(stateenergy_elastic[q_point]) * sqrt(abs(stateenergy_elastic[q_point]) / mass)
                                ΔV₂ += -sign(V) * sign(stateenergy_plastic[q_point]) * sqrt(abs(stateenergy_plastic[q_point]) / mass)
                                # ΔV₂ += -sign(V) * sqrt(stateenergy_elastic[q_point] * abs(getnormal(facetvalues, q_point) ⋅ (sign(V) .* [0, 0, 1])) / mass)
                            end
                            ΔV₁ += ΔV + ΔV₂ / 4
                            # ΔV₁ += ΔV + ΔV₂
                        end
                        ΔV₁ /= length(getcellset(grid, "end_plate"))
                        ΔV = ΔV₁
                        V = -V
                    end
                    state_strain_itr       .= state_strain
                    state_energy_total_itr .= state_energy_total
                end # end of itx
            end


            # break if within tolerance
            norm_r = norm(r[Ferrite.free_dofs(ch)])
            @printf("\tIteration: %02d \tresidual: %.9f\n", newton_itr, norm_r) # \titx coords: $intersecting_coords") # , $(@sprintf("%.9f", norm_s))")
            if norm_r < NEWTON_TOL
                break
            end


            # increment residual to next NR iteration
            apply_zero!(K, r, ch)
            Δu = Symmetric(K) \ r
            apply_zero!(Δu, ch)
            u -= Δu
            newton_itr += 1
        end
        if newton_itr > NEWTON_M
            error("Reached maximum Newton iterations, aborting")
        end

        if do_verify
            # ϵ̇ .= (u - u₀) / Δt # (dN * Δt)
            # if use_rve_or_E8
            #     i = 8; println_i_rve(i, state_traction, state_energy_total, state_energy_elastic, state_energy_plastic, p, ΔV)
            # else
            #     i = 23; println_i_E8(i, state_traction, state_energy_total, state_energy_elastic, state_energy_plastic, p, ΔV)
            # end
            # # println((get_node_coordinate(gridnodes[i]), u[vectorstrafe_i(i)]))
        else
            # if 80 < timestep < 140
            # end
            i = node_sphere_o; println_i(i, state_traction, state_energy_total, state_energy_elastic, state_energy_plastic, p, ΔV)
        end


        # write to results file
            # velocity = zeros(ncells_grid)
            strainᵗ_vonMises    = zeros(gridncells)
            strainᵗ_11          = zeros(gridncells)
            strainᵗ_22          = zeros(gridncells)
            strainᵗ_33          = zeros(gridncells)
            strainᵗ_12          = zeros(gridncells)
            strainᵗ_13          = zeros(gridncells)
            strainᵗ_21          = zeros(gridncells)
            strainᵖ_vonMises    = zeros(gridncells)
            strainᵖ_11          = zeros(gridncells)
            strainᵖ_22          = zeros(gridncells)
            strainᵖ_33          = zeros(gridncells)
            strainᵖ_12          = zeros(gridncells)
            strainᵖ_13          = zeros(gridncells)
            strainᵖ_21          = zeros(gridncells)
            strainᵉ_vonMises    = zeros(gridncells)
            strainᵉ_11          = zeros(gridncells)
            strainᵉ_22          = zeros(gridncells)
            strainᵉ_33          = zeros(gridncells)
            strainᵉ_12          = zeros(gridncells)
            strainᵉ_13          = zeros(gridncells)
            strainᵉ_21          = zeros(gridncells)
            strainrateᵗ_vonMises= zeros(gridncells)
            strainrateᵗ_11      = zeros(gridncells)
            strainrateᵗ_22      = zeros(gridncells)
            strainrateᵗ_33      = zeros(gridncells)
            strainrateᵗ_12      = zeros(gridncells)
            strainrateᵗ_13      = zeros(gridncells)
            strainrateᵗ_21      = zeros(gridncells)
            strainrateᵖ_vonMises= zeros(gridncells)
            strainrateᵖ_11      = zeros(gridncells)
            strainrateᵖ_22      = zeros(gridncells)
            strainrateᵖ_33      = zeros(gridncells)
            strainrateᵖ_12      = zeros(gridncells)
            strainrateᵖ_13      = zeros(gridncells)
            strainrateᵖ_21      = zeros(gridncells)
            strainrateᵉ_vonMises= zeros(gridncells)
            strainrateᵉ_11      = zeros(gridncells)
            strainrateᵉ_22      = zeros(gridncells)
            strainrateᵉ_33      = zeros(gridncells)
            strainrateᵉ_12      = zeros(gridncells)
            strainrateᵉ_13      = zeros(gridncells)
            strainrateᵉ_21      = zeros(gridncells)
            sigma_vonMises      = zeros(gridncells)
            sigma_11            = zeros(gridncells)
            sigma_22            = zeros(gridncells)
            sigma_33            = zeros(gridncells)
            sigma_12            = zeros(gridncells)
            sigma_13            = zeros(gridncells)
            sigma_21            = zeros(gridncells)
            # κ_values = zeros(gridncells)
            # for (el, state_cells) in enumerate(eachcol(state_material))
            for ((el, state_cells), state_old_cells) in zip(enumerate(eachcol(state_material)), eachcol(state_material_old))
                for (state, state_old) in zip(state_cells, state_old_cells)
                    # velocity[el] += first(rand(n_timesteps))
                    strainᵗ_vonMises[el]    += vonMises(state.ϵ)
                    strainᵗ_11[el]          += state.ϵ[1, 1]
                    strainᵗ_22[el]          += state.ϵ[2, 2]
                    strainᵗ_33[el]          += state.ϵ[3, 3]
                    strainᵗ_12[el]          += state.ϵ[1, 2]
                    strainᵗ_13[el]          += state.ϵ[1, 3]
                    strainᵗ_21[el]          += state.ϵ[2, 1]
                    strainᵖ_vonMises[el]    += vonMises(state.ϵᵖ)
                    strainᵖ_11[el]          += state.ϵᵖ[1, 1]
                    strainᵖ_22[el]          += state.ϵᵖ[2, 2]
                    strainᵖ_33[el]          += state.ϵᵖ[3, 3]
                    strainᵖ_12[el]          += state.ϵᵖ[1, 2]
                    strainᵖ_13[el]          += state.ϵᵖ[1, 3]
                    strainᵖ_21[el]          += state.ϵᵖ[2, 1]
                    strainᵉ_vonMises[el]    += (vonMises(state.ϵ) - vonMises(state.ϵᵖ))
                    strainᵉ_11[el]          += (state.ϵ[1, 1] - state.ϵᵖ[1, 1])
                    strainᵉ_22[el]          += (state.ϵ[2, 2] - state.ϵᵖ[2, 2])
                    strainᵉ_33[el]          += (state.ϵ[3, 3] - state.ϵᵖ[3, 3])
                    strainᵉ_12[el]          += (state.ϵ[1, 2] - state.ϵᵖ[1, 2])
                    strainᵉ_13[el]          += (state.ϵ[1, 3] - state.ϵᵖ[1, 3])
                    strainᵉ_21[el]          += (state.ϵ[2, 1] - state.ϵᵖ[2, 1])
                    strainrateᵗ_vonMises[el]+= (vonMises(state.ϵ) - vonMises(state_old.ϵ)) / Δt # (dN * Δt)
                    strainrateᵗ_11[el]      += (state.ϵ[1, 1] - state_old.ϵ[1, 1]) / Δt # (dN * Δt)
                    strainrateᵗ_22[el]      += (state.ϵ[2, 2] - state_old.ϵ[2, 2]) / Δt # (dN * Δt)
                    strainrateᵗ_33[el]      += (state.ϵ[3, 3] - state_old.ϵ[3, 3]) / Δt # (dN * Δt)
                    strainrateᵗ_12[el]      += (state.ϵ[1, 2] - state_old.ϵ[1, 2]) / Δt # (dN * Δt)
                    strainrateᵗ_13[el]      += (state.ϵ[1, 3] - state_old.ϵ[1, 3]) / Δt # (dN * Δt)
                    strainrateᵗ_21[el]      += (state.ϵ[2, 1] - state_old.ϵ[2, 1]) / Δt # (dN * Δt)
                    strainrateᵖ_vonMises[el]+= (vonMises(state.ϵᵖ) - vonMises(state_old.ϵᵖ)) / Δt # (dN * Δt)
                    strainrateᵖ_11[el]      += (state.ϵᵖ[1, 1] - state_old.ϵᵖ[1, 1]) / Δt # (dN * Δt)
                    strainrateᵖ_22[el]      += (state.ϵᵖ[2, 2] - state_old.ϵᵖ[2, 2]) / Δt # (dN * Δt)
                    strainrateᵖ_33[el]      += (state.ϵᵖ[3, 3] - state_old.ϵᵖ[3, 3]) / Δt # (dN * Δt)
                    strainrateᵖ_12[el]      += (state.ϵᵖ[1, 2] - state_old.ϵᵖ[1, 2]) / Δt # (dN * Δt)
                    strainrateᵖ_13[el]      += (state.ϵᵖ[1, 3] - state_old.ϵᵖ[1, 3]) / Δt # (dN * Δt)
                    strainrateᵖ_21[el]      += (state.ϵᵖ[2, 1] - state_old.ϵᵖ[2, 1]) / Δt # (dN * Δt)
                    strainrateᵉ_vonMises[el]+= ((vonMises(state.ϵ) - vonMises(state.ϵᵖ)) - (vonMises(state_old.ϵ) - vonMises(state_old.ϵᵖ))) / Δt # (dN * Δt)
                    strainrateᵉ_11[el]      += ((state.ϵ[1, 1] - state.ϵᵖ[1, 1]) - (state_old.ϵ[1, 1] - state_old.ϵᵖ[1, 1])) / Δt # (dN * Δt)
                    strainrateᵉ_22[el]      += ((state.ϵ[2, 2] - state.ϵᵖ[2, 2]) - (state_old.ϵ[2, 2] - state_old.ϵᵖ[2, 2])) / Δt # (dN * Δt)
                    strainrateᵉ_33[el]      += ((state.ϵ[3, 3] - state.ϵᵖ[3, 3]) - (state_old.ϵ[3, 3] - state_old.ϵᵖ[3, 3])) / Δt # (dN * Δt)
                    strainrateᵉ_12[el]      += ((state.ϵ[1, 2] - state.ϵᵖ[1, 2]) - (state_old.ϵ[1, 2] - state_old.ϵᵖ[1, 2])) / Δt # (dN * Δt)
                    strainrateᵉ_13[el]      += ((state.ϵ[1, 3] - state.ϵᵖ[1, 3]) - (state_old.ϵ[1, 3] - state_old.ϵᵖ[1, 3])) / Δt # (dN * Δt)
                    strainrateᵉ_21[el]      += ((state.ϵ[2, 1] - state.ϵᵖ[2, 1]) - (state_old.ϵ[2, 1] - state_old.ϵᵖ[2, 1])) / Δt # (dN * Δt)
                    sigma_vonMises[el]      += vonMises(state.σ)
                    sigma_11[el]            += state.σ[1, 1]
                    sigma_22[el]            += state.σ[2, 2]
                    sigma_33[el]            += state.σ[3, 3]
                    sigma_12[el]            += state.σ[1, 2]
                    sigma_13[el]            += state.σ[1, 3]
                    sigma_21[el]            += state.σ[2, 1]
                    # κ_values[el] += state.k*material_endplate.H
                end
                # velocity[el] /= length(cell_states)
                strainᵗ_vonMises[el]    /= length(state_cells)
                strainᵗ_11[el]          /= length(state_cells)
                strainᵗ_22[el]          /= length(state_cells)
                strainᵗ_33[el]          /= length(state_cells)
                strainᵗ_12[el]          /= length(state_cells)
                strainᵗ_13[el]          /= length(state_cells)
                strainᵗ_21[el]          /= length(state_cells)
                strainᵖ_vonMises[el]    /= length(state_cells)
                strainᵖ_11[el]          /= length(state_cells)
                strainᵖ_22[el]          /= length(state_cells)
                strainᵖ_33[el]          /= length(state_cells)
                strainᵖ_12[el]          /= length(state_cells)
                strainᵖ_13[el]          /= length(state_cells)
                strainᵖ_21[el]          /= length(state_cells)
                strainᵉ_vonMises[el]    /= length(state_cells)
                strainᵉ_11[el]          /= length(state_cells)
                strainᵉ_22[el]          /= length(state_cells)
                strainᵉ_33[el]          /= length(state_cells)
                strainᵉ_12[el]          /= length(state_cells)
                strainᵉ_13[el]          /= length(state_cells)
                strainᵉ_21[el]          /= length(state_cells)
                strainrateᵗ_vonMises[el]/= length(state_cells)
                strainrateᵗ_11[el]      /= length(state_cells)
                strainrateᵗ_22[el]      /= length(state_cells)
                strainrateᵗ_33[el]      /= length(state_cells)
                strainrateᵗ_12[el]      /= length(state_cells)
                strainrateᵗ_13[el]      /= length(state_cells)
                strainrateᵗ_21[el]      /= length(state_cells)
                strainrateᵖ_vonMises[el]/= length(state_cells)
                strainrateᵖ_11[el]      /= length(state_cells)
                strainrateᵖ_22[el]      /= length(state_cells)
                strainrateᵖ_33[el]      /= length(state_cells)
                strainrateᵖ_12[el]      /= length(state_cells)
                strainrateᵖ_13[el]      /= length(state_cells)
                strainrateᵖ_21[el]      /= length(state_cells)
                strainrateᵉ_vonMises[el]/= length(state_cells)
                strainrateᵉ_11[el]      /= length(state_cells)
                strainrateᵉ_22[el]      /= length(state_cells)
                strainrateᵉ_33[el]      /= length(state_cells)
                strainrateᵉ_12[el]      /= length(state_cells)
                strainrateᵉ_13[el]      /= length(state_cells)
                strainrateᵉ_21[el]      /= length(state_cells)
                sigma_vonMises[el]      /= length(state_cells)
                sigma_11[el]            /= length(state_cells)
                sigma_22[el]            /= length(state_cells)
                sigma_33[el]            /= length(state_cells)
                sigma_12[el]            /= length(state_cells)
                sigma_13[el]            /= length(state_cells)
                sigma_21[el]            /= length(state_cells)
                # κ_values[el] /= length(state_cells)
            end
            strain_max[timestep] = maximum(abs, strainᵗ_vonMises) # maximum displacement in current timestep
            strainrate_max[timestep] = if do_verify && !use_rve_or_E8
                try
                    u_gaugeright = average([u[vectorstrafe_i(n)][1] for n in getnodeset(grid, "gauge_right")])
                    u_gaugeleft = average([u[vectorstrafe_i(n)][1] for n in getnodeset(grid, "gauge_left")])
                    (u_gaugeright - u_gaugeleft) / Δt
                catch
                    coord_gaugeright = average([get_node_coordinate(grid, n)[1] for n in getnodeset(grid, "gauge_right")])
                    coord_gaugeleft = average([get_node_coordinate(grid, n)[1] for n in getnodeset(grid, "gauge_left")])
                    ((coord_gaugeright - coord_gaugeleft) / gauge_length) / Δt
                end
            else
                maximum(abs, strainrateᵗ_vonMises) # maximum displacement in current timestep
            end
            stress_max[timestep] = maximum(abs, sigma_vonMises) # maximum displacement in current timestep
            # strain_max[timestep] = maximum(abs, strainᵗ_22) # maximum displacement in current timestep
            # stress_max[timestep] = maximum(abs, sigma_22) # maximum displacement in current timestep
            Vₐ = try
                derive_serial(collect(range(0, 0.1, n_timesteps))[begin:Int64((0.02 / (0.1 / (n_timesteps - 1))) ÷ 1 + 1)][begin:timestep],
                    ΔX, 1:timestep, timestep, :five, 1, Δt)
            catch
                try
                    (ΔX[timestep] - ΔX[timestep - 1]) / Δt
                catch
                    V
                end
            end
            velocity_applied = zeros((3, gridnnodes))
            velocity_actual = zeros((3, gridnnodes))
            if do_verify
                for i in getnodeset(grid, use_rve_or_E8 ? "front_xz" : "right")
                    velocity_applied[:, i] = V .* [use_rve_or_E8 ? 1.0 : 0.0, use_rve_or_E8 ? 1.0 : 0.0, 0.0]
                    velocity_actual[:, i] = Vₐ .* [use_rve_or_E8 ? 1.0 : 0.0, use_rve_or_E8 ? 1.0 : 0.0, 0.0]
                end
            else
                for i in getnodeset(grid, "indenter_tip")
                    velocity_applied[:, i] = V .* [0., 0., 1.]
                    velocity_actual[:, i] = Vₐ .* [0., 0., 1.]
                end
            end
            # velocity = zeros(n_dofs)
            # for i in getnodeset(grid, "indenter_tip")
            #     velocity[vectorstrafe_i(i)] .= V .* [0., 0., 1.]
            # end
            VTKGridFile(filename * "-t$timestep", dh) do vtk
                # vtk["TimeValue"] = float(timestep)
                write_solution(vtk, dh, u) # displacement field
                # // DONE [20250110] (JMA3): figure out how to edit output vectors for visualization
                # [20250206T0911] (JMA3): `write_node_data` ≡ `write_solution` for display vectors in ParaView
                # write_solution(vtk, dh, fill(t, n_dofs), "v") # velocity field
                # write_cell_data(vtk, velocity, "v") # velocity field
                write_node_data(vtk, velocity_applied, "Velocity (Applied)")
                write_node_data(vtk, velocity_actual, "Velocity (Actual)")
                # write_solution(vtk, dh, u) # displacement field

                write_cell_data(vtk, strainᵗ_vonMises,      "strain_t")
                write_cell_data(vtk, strainᵗ_11,            "strain_t_11")
                write_cell_data(vtk, strainᵗ_22,            "strain_t_22")
                write_cell_data(vtk, strainᵗ_33,            "strain_t_33")
                write_cell_data(vtk, strainᵗ_12,            "strain_t_12")
                write_cell_data(vtk, strainᵗ_13,            "strain_t_13")
                write_cell_data(vtk, strainᵗ_21,            "strain_t_21")
                write_cell_data(vtk, strainᵖ_vonMises,      "strain_p")
                write_cell_data(vtk, strainᵖ_11,            "strain_p_11")
                write_cell_data(vtk, strainᵖ_22,            "strain_p_22")
                write_cell_data(vtk, strainᵖ_33,            "strain_p_33")
                write_cell_data(vtk, strainᵖ_12,            "strain_p_12")
                write_cell_data(vtk, strainᵖ_13,            "strain_p_13")
                write_cell_data(vtk, strainᵖ_21,            "strain_p_21")
                write_cell_data(vtk, strainᵉ_vonMises,      "strain_e")
                write_cell_data(vtk, strainᵉ_11,            "strain_e_11")
                write_cell_data(vtk, strainᵉ_22,            "strain_e_22")
                write_cell_data(vtk, strainᵉ_33,            "strain_e_33")
                write_cell_data(vtk, strainᵉ_12,            "strain_e_12")
                write_cell_data(vtk, strainᵉ_13,            "strain_e_13")
                write_cell_data(vtk, strainᵉ_21,            "strain_e_21")
                write_cell_data(vtk, strainrateᵗ_vonMises,  "strainrate_t")
                write_cell_data(vtk, strainrateᵗ_11,        "strainrate_t_11")
                write_cell_data(vtk, strainrateᵗ_22,        "strainrate_t_22")
                write_cell_data(vtk, strainrateᵗ_33,        "strainrate_t_33")
                write_cell_data(vtk, strainrateᵗ_12,        "strainrate_t_12")
                write_cell_data(vtk, strainrateᵗ_13,        "strainrate_t_13")
                write_cell_data(vtk, strainrateᵗ_21,        "strainrate_t_21")
                write_cell_data(vtk, strainrateᵖ_vonMises,  "strainrate_p")
                write_cell_data(vtk, strainrateᵖ_11,        "strainrate_p_11")
                write_cell_data(vtk, strainrateᵖ_22,        "strainrate_p_22")
                write_cell_data(vtk, strainrateᵖ_33,        "strainrate_p_33")
                write_cell_data(vtk, strainrateᵖ_12,        "strainrate_p_12")
                write_cell_data(vtk, strainrateᵖ_13,        "strainrate_p_13")
                write_cell_data(vtk, strainrateᵖ_21,        "strainrate_p_21")
                write_cell_data(vtk, strainrateᵉ_vonMises,  "strainrate_e")
                write_cell_data(vtk, strainrateᵉ_11,        "strainrate_e_11")
                write_cell_data(vtk, strainrateᵉ_22,        "strainrate_e_22")
                write_cell_data(vtk, strainrateᵉ_33,        "strainrate_e_33")
                write_cell_data(vtk, strainrateᵉ_12,        "strainrate_e_12")
                write_cell_data(vtk, strainrateᵉ_13,        "strainrate_e_13")
                write_cell_data(vtk, strainrateᵉ_21,        "strainrate_e_21")
                write_cell_data(vtk, sigma_vonMises,        "sigma_t")
                write_cell_data(vtk, sigma_11,              "sigma_11")
                write_cell_data(vtk, sigma_22,              "sigma_22")
                write_cell_data(vtk, sigma_33,              "sigma_33")
                write_cell_data(vtk, sigma_12,              "sigma_12")
                write_cell_data(vtk, sigma_13,              "sigma_13")
                write_cell_data(vtk, sigma_21,              "sigma_21")
                # for i in 1:3, j in 1:3
                #     σij = [x[i, j] for x in σ]
                #     write_cell_data(vtk, σij, "sigma_$(i)$(j)")
                # end
                # write_cell_data(vtk, κ_values, "Drag stress [Pa]")
                # Ferrite.write_constraints(vtk, ch)
                # pvd[timestep] = vtk
                collection_add_timestep(pvd, vtk, t)
            end # closes vtk
        #


        # update the old states with the converged values for next timestep
        u₀ = u
        V += ΔV
        # states_material_old        .= state_material
        # states_strain_rate         .= (state_strain - state_strain_old) / Δt # (dN * Δt)
        # states_strain_old          .= state_strain
        # # states_energy_total_old    .= state_energy_total
        # # states_energy_elastic_old  .= state_energy_elastic
        # # states_energy_plastic_old  .= state_energy_plastic
        # states_energy_total[timestep]   = state_energy_total
        # states_energy_elastic[timestep] = state_energy_elastic
        # states_energy_plastic[timestep] = state_energy_plastic
        states_material[timestep]       = state_material
        states_strain_rate[timestep]    = (state_strain - state_strain_old) / Δt # (dN * Δt)
        states_energy_total[timestep]   = state_energy_total
        states_energy_elastic[timestep] = state_energy_elastic
        states_energy_plastic[timestep] = state_energy_plastic
        u_max[timestep] = maximum(abs, u) # maximum displacement in current timestep
    end


    # save results and exit
    println("\nTime to solve: ", timeduration(simulationtime_0, Dates.now()))
    vtk_save(pvd);
    # return (u_max, v_max), traction_magnitude, (u_max ./ plate_thickness, v_max ./ R), stress_max
    return do_verify ? (strain_max, strainrate_max, stress_max) : nothing
end

if do_verify
    strain, strainrate, stress = solve();
    # strain, strainrate, stress = strain[begin:2:end-1], stress[begin:2:end-1]
    println("Size of vectors to plot: ", map(size, [expdata_strain, expdata_stress, modeldata_strain, modeldata_stress, strain, stress]))
    rmse_model = (!do_verify || !use_rve_or_E8) ? NaN : rmse((expdata_strain, expdata_stress), (modeldata_strain, modeldata_stress))
    rmse_verify = (!do_verify || !use_rve_or_E8) ? NaN : rmse((expdata_strain, expdata_stress), (strain, stress))
    default(fontfamily="Computer Modern", grid=false, framestyle=:box)
    p = plot(
        [expdata_strain, modeldata_strain, strain],
        [expdata_stress, modeldata_stress, stress],
        label=["Raw" @sprintf("Point, E: %.03f [MPa]", rmse_model) @sprintf("Ferrite.jl (%s), E: %.03f [MPa]", (use_rve_or_E8 ? "RVE" : "E8"), rmse_verify)],
        markershape=:auto,
        linewidth=2,
        legendposition=:bottomright,
        xlabel="Strain [mm/mm]",
        ylabel="Stress [MPa]",
    )
    percentdifference_model_verify = average(map(percentdifference, stress[map(x->(y = findfirst(x .<= strain); !isnothing(y) ? y : findlast(x .>= strain)), modeldata_strain)], modeldata_stress))
    annotate!(p, 0.15, 100.0, text(@sprintf("%3.3f %%\$_{diff}\$ between\nPoint and Ferrite.jl", percentdifference_model_verify), "Computer Modern", :left, 8))
    display(p)
    # rmse_model = (!do_verify || !use_rve_or_E8) ? NaN : rmse((expdata_strain, expdata_strainrate), (modeldata_strain, modeldata_strainrate))
    # rmse_verify = (!do_verify || !use_rve_or_E8) ? NaN : rmse((expdata_strain, expdata_strainrate), (strain, strainrate))
    # q = plot(
    #     [expdata_strain, modeldata_strain, strain],
    #     [expdata_strainrate, modeldata_strainrate, strainrate],
    #     label=["Raw" @sprintf("Point, E: %3.3e [s\$^{-1}\$]", rmse_model) @sprintf("Ferrite.jl (%s), E: %3.3e [s\$^{-1}\$]", (use_rve_or_E8 ? "RVE" : "E8"), rmse_verify)],
    #     markershape=:auto,
    #     linewidth=2,
    #     xlabel="Strain [mm/mm]",
    #     ymirror=true,
    #     ylabel="Strain-Rate [mm/mm/s]",
    # )
    # r = plot(p, q, layout=(1, 2))
    # savefig(r, "_plot.png")
    # display(r)
else
    solve();
end