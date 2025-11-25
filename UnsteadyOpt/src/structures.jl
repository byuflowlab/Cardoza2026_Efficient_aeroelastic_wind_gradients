


"""
    get_constant_distributedload(element, fy, fz)
Assume a constant distributed load along the length of the element, integrate and prepare to apply to the points. 
This assumes no moment and no load in the axial direction. 

*Arguments*
- `element`: The element to which the load is applied.
- `fy`: The distributed load in the y direction.
- `fz`: The distributed load in the z direction.

*Returns*
- `GXBeam.DistributedLoads`: The distributed load object with the specified loads.
"""
function get_constant_distributedload(element, fy, fz)
    L = element.L
    TF = promote_type(typeof(fy), typeof(fz))  
    zero_load = StaticArrays.SVector{3, TF}(0.0, 0.0, 0.0)
    return GXBeam.DistributedLoads{TF}(StaticArrays.SVector{3, TF}(0.0, fy*L/2, fz*L/2), StaticArrays.SVector{3, TF}(0.0, fy*L/2, fz*L/2), zero_load, zero_load, zero_load, zero_load, zero_load, zero_load)
end


function compute_blade_mass(mass_list, elements, TF)
    total_mass = TF(0.0)
    for i in eachindex(elements)
        element_mass = mass_list[i][1] * elements[i].L
        total_mass += element_mass
    end
    return total_mass
end




function get_steady_parameters(p, points, start, stop, xp, n_compliance, n_mass, nelem)
    # compliance matrix for each beam element
    compliance = create_list_property(p[1:n_compliance])
    mass = create_list_property(p[n_compliance+1:n_compliance+n_mass])
    loads = p[n_compliance+n_mass+1:end]
    fy_i = loads[1:nelem]  
    fz_i = loads[nelem+1:end]

    # create assembly of interconnected nonlinear beams
    assembly = GXBeam.Assembly(points, start, stop, compliance=compliance, midpoints=xp, mass=mass)
    
    #Update the distributed loads.  
    distributed_loads = Dict(i => get_constant_distributedload(assembly.elements[i], fy_i[i], fz_i[i]) for i in 1:nelem)
    
    # return named tuple with new arguments
    return (; assembly, distributed_loads)
end


function extreme_loading_analysis!(g, clt_list, compliance, mass, fy, fz, points, xp, azimuth, eps_ult, gravity_value, num_elements, num_buckling, N_buckling, N_elements, Nshift, buckling_scale, strain_scale, deflection_scale)

    ### Store information in the p vector for the adjoint
    p_compliance = get_list_property(compliance)
    p_mass = get_list_property(mass)
    p_loads = vcat(fy, fz)
    p = vcat(p_compliance, p_mass, p_loads)

    # @show typeof(p)

    n_compliance = length(p_compliance)
    n_mass = length(p_mass)
    

    np = length(points)
    nelem = length(fy) 

    pfunc = (pee, t) -> get_steady_parameters(
        pee,
        points,
        1:nelem,
        2:np,
        xp,
        n_compliance,
        n_mass,
        nelem) 
        
    # if isa(p[1], ReverseDiff.TrackedReal) #TODO: Probably ought to put this in a function of some sort. 
    #     assembly = pfunc(ReverseDiff.value.(p), 0.0).assembly
    # else
    #     assembly = pfunc(ForwardDiff.value.(p), 0.0).assembly #Get the assembly from the pfunc.
    # end
    assembly = pfunc(nondual_value.(p), 0.0).assembly #Get the assembly from the pfunc.


    gravity = SVector(-gravity_value*cos(azimuth), -gravity_value*sin(azimuth), 0.0)
     
    prescribed_conditions = Dict(1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Fixed first node. 

    # error("Stop here. ")
    ### Run static analysis
    _, state, converged = GXBeam.steady_state_analysis(assembly; prescribed_conditions, pfunc, p, gravity, iterations=100, ftol=1e-8) 
    
    #A flag to print something if the analysis didn't converge.
    if !converged
        println("Analysis did not converge.")
    end
    

    ### Fill constraints  
    for i in 1:nelem
        # Fi, Mi = rotate_internal_loads(state.elements[i].Fi, state.elements[i].Mi, 0.0)
        Fi = state.elements[i].Fi
        Mi = state.elements[i].Mi
        
        fb, fs, _ = check_failure(Fi, Mi, clt_list[i], eps_ult)

        b_idx, fs_idx = get_failure_constraint_indices(i, num_buckling, num_elements, N_buckling)

        b_idx = b_idx .+ Nshift
        fs_idx = fs_idx .+ Nshift
        # @show b_idx, fs_idx

        g[b_idx] = fb./buckling_scale #Buckling constraint
        g[fs_idx] = fs/strain_scale #Strain constraint
    end  

    ### Deflection constraint
    g[Nshift + N_buckling + N_elements + 1] = state.elements[end].u[3]/deflection_scale
    # @show Nshift + N_buckling + N_elements + 1
end



function get_unsteady_parameters(p, points, start, stop, xp, azimuth0, g, n_compliance, n_mass, nelem)
    
    compliance = create_list_property(p[1:n_compliance])
    mass = create_list_property(p[n_compliance+1:n_compliance+n_mass])
    loads = p[n_compliance+n_mass+1:n_compliance+n_mass+(2*nelem)]
    fy_i = loads[1:nelem]  
    fz_i = loads[nelem+1:end]
    omega_internal = p[n_compliance+n_mass+(2*nelem)+1] 

    # create assembly of interconnected nonlinear beams
    assembly = GXBeam.Assembly(points, start, stop, compliance=compliance, midpoints=xp, mass=mass)

    #Update the distributed loads.  
    distributed_loads = Dict(i => get_constant_distributedload(assembly.elements[i], fy_i[i], fz_i[i]) for i in 1:nelem)

    gravity = (tee) -> SVector(-g*cos(azimuth0 + omega_internal*tee), -g*sin(azimuth0 + omega_internal*tee), 0.0)

    angular_velocity = SVector(0.0, 0.0, omega_internal)
    
    # return named tuple with new arguments
    return (; assembly, distributed_loads, gravity, angular_velocity)
end



function fatigue_analysis!(g, clt_list, compliance, mass, chords, twists, Omega_rated, rated_pitch, points, xp, tvec, B, rvec, rhub, rtip, xcp, airfoils, hubHt, azimuth0, yaw, tilt, precone, Ufit, Vfit, Wfit, shearExp, rho, mu, a, eps_ult, m, nu, gravity_value, fat_idxs, num_elements, N_bulk_full, Nshift, Ntimecon, deflection_scale)

    TF = typeof(chords[1])

    p_compliance = get_list_property(compliance)
    p_mass = get_list_property(mass)
    
    if isa(chords[1], ReverseDiff.TrackedReal)
        # @show typeof(chords)
        p = get_vector_tracked_reals(p_compliance, p_mass, chords, twists, Omega_rated, rated_pitch) #Hopefully this doesn't need to be a TrackedArray. 
        # p = get_tracked_array(chords, p_compliance, p_mass, chords, twists, Omega_rated, rated_pitch) #The first argurment is just to get the type.
        # @show typeof(p)
    else
        p = vcat(p_compliance, p_mass, chords, twists, Omega_rated, rated_pitch) #Original line.
    end

    n_compliance = length(p_compliance)
    n_mass = length(p_mass)

    np = length(points)
    nelem = length(chords)

    pfunc = (pee, tee) -> get_unsteady_parameters(
        pee,
        points,
        1:nelem,
        2:np,
        xp,
        azimuth0,
        gravity_value,
        n_compliance,
        n_mass,
        nelem)

    if isa(chords[1], ReverseDiff.TrackedReal)
        # @show typeof(p), typeof(p[1])
        # @show ReverseDiff.value(p[1])
        # @show ReverseDiff.value.(p)
        assembly = pfunc(ReverseDiff.value.(p), 0.0).assembly #Get the assembly from the pfunc.
    else
         assembly = pfunc(ForwardDiff.value.(p), 0.0).assembly #Get the assembly from the pfunc.
    end
    # assembly = pfunc(ForwardDiff.value.(p), 0.0).assembly #Get the assembly from the pfunc.
    # @show typeof(assembly)

    ufun(t) = SVector(Ufit(t), Vfit(t), Wfit(t))
    omegafun(t) = SVector(0.0, 0.0, 0.0)
    udotfun(t) = SVector(0.0, 0.0, 0.0)
    omegadotfun(t) = SVector(0.0, 0.0, 0.0)
    Vinf(t) = Ufit(t)
    RS(t) = Omega_rated
    Vinfdot(t) = 0.0
    RSdot(t) = 0.0
    env = WATT.SimpleEnvironment(rho, mu, a, shearExp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)

    rotor = WATT.Rotor(Int(B), hubHt, true; tilt, yaw)
    # @show typeof(chords)
    # @show typeof(twists)
    blade = WATT.Blade(rvec, chords, twists, xcp, airfoils; rhub=rhub, rtip=rtip, precone)
    
    aerostates, gxhistory, mesh = WATT.initialize_sim(blade, assembly, tvec; verbose=false, pfunc=pfunc, p=p)

    function prepp!(p, Fx, Fy, mx) #Todo: I'm not sure that I have the load transforms correct here. What should be fed in here? should it be rotated here or is this assuming the rotated loads? 
        loads = @views p[n_compliance+n_mass+1:n_compliance+n_mass+(2*nelem)]
        loads[1:nelem] .= Fy #fy
        loads[nelem+1:end] .= -Fx #fz
    end

    WATT.run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false, prepp=prepp!, p=p)


    ntime = length(tvec)

    ### Extract the tip deflection
    tipdef_idxs = Nshift+1:Nshift+Ntimecon
    for i = 201:ntime
        g[tipdef_idxs[i-200]] = gxhistory[i].elements[end].u[3]/deflection_scale
    end
    # @show tipdef_idxs



    ### Extract loads at the locations we want to check for damage.
    forces = zeros(TF, length(fat_idxs), ntime, 3) #todo: instead of allocating this... I could probably just pass the gxhistory directly to the damage calculation function and extract the forces and moments there.
    moments = zeros(TF, length(fat_idxs), ntime, 3)
    for i in eachindex(tvec)
        for j in eachindex(fat_idxs)
            forces[j, i, :] = gxhistory[i].elements[fat_idxs[j]].Fi
            moments[j, i, :] = gxhistory[i].elements[fat_idxs[j]].Mi
        end
    end


    t_elapsed = tvec[end]-tvec[1]
    
    num_constraints = Nshift+Ntimecon #Number of constraints up until this point in the optimization. 

    calculate_damage(
        g,
        fat_idxs,
        forces,
        moments,
        clt_list,
        eps_ult,
        m,
        t_elapsed,
        nu,
        num_constraints,
        Omega_rated,
        num_elements)
end

