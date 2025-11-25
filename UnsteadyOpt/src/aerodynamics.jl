


function allocate_interp_airfoil(i, path)
    M = readdlm(joinpath(path, "airfoil_$i.dat"), ',')
    # Note: Inherently type unstable. 
    if i <= 13
        return CCBlade.Cylinder(maximum(M[:,3]))
    else
        return CCBlade.AlphaAF(M[:,1], M[:,2], M[:,3])
    end
end

function get_interp_polars(rvec, path)
    return [allocate_interp_airfoil(i, path) for i in eachindex(rvec)]
end


function get_loads(outs, chords, rho, TF)
    nr = length(outs)

    fy = zeros(TF, nr)
    fz = zeros(TF, nr)

    for i = 1:nr
        sphi, cphi = sincos(outs[i].phi)
        Cx = outs[i].cl*cphi + outs[i].cd*sphi #Rotate into the blade root reference frame. I just checked this. 6/11/25
        Cy = -(outs[i].cl*sphi - outs[i].cd*cphi)
        q_local = 0.5*rho*outs[i].W^2 #Local dynamic pressure

        #Dimensionalize the loads and rotate into GXBeam root frame. 
        fy[i] = Cy*q_local*chords[i]
        fz[i] = -Cx*q_local*chords[i]
    end

    return fy, fz
end



function get_thrustpowercurve(rotor, sections, Vcurve, tsr, pitches, yaw, tilt, azimuth, hubHt, shearExp, rho, Vmean, TF)

    precone = rotor.precone

    maxomega = 12.0*2*pi/60
    rotorR = rotor.Rtip*cos(precone)
    rvec = sections.r 

    powercurve = zeros(TF, length(Vcurve))
    Tcurve = zeros(TF, length(Vcurve))
    

    for i in eachindex(Vcurve)
        omega_i = min(Vcurve[i]*tsr/rotorR, maxomega) #Calculating rotational speed and constraining it to a maximum value.
        pitch_i = pitches[i] 
        opi = windturbine_op.(Vcurve[i], omega_i, pitch_i, rvec, precone, yaw, tilt, azimuth, hubHt, shearExp, rho) 
        out_i = CCBlade.solve.(Ref(rotor), sections, opi) 
        T_i, Q_i = thrusttorque(rotor, sections, out_i) 
        powercurve[i] = Q_i*omega_i
        Tcurve[i] = T_i
    end

    kay = 2.0  # weibull shape
    A = Vmean / gamma(1.0 + 1.0/kay)
    cdf = 1.0 .- exp.(-(Vcurve/A).^kay) 
    AEP = FLOWMath.trapz(cdf, powercurve) * 365*24

    return Tcurve, powercurve, AEP
end


