#=
Time the reverse mode AD while varying the number of sparse design variables.

Adam Cardoza 10/2/25
=# 
localpath = @__DIR__
cd(localpath)

using GXBeamCS, GXBeam, CCBlade, OpenFASTTools, DynamicStallModels
using FLOWMath, DelimitedFiles, LinearAlgebra, Dates
using StaticArrays, StructArrays, SparseArrays
using ForwardDiff, FiniteDiff, DiffResults, PolyesterForwardDiff, ReverseDiff
using SNOW, Snopt
using Plots, LaTeXStrings
using UnsteadyOpt

df = Dates.DateFormat("yymmdd_HH.MM.SS")
now = Dates.now()
nowstr = Dates.format(now, df)
filename = splitpath(@__FILE__)[end]
rootname = "_checkmemory_"


println("running ", filename, " at ", nowstr)

of = OpenFASTTools
uo = UnsteadyOpt
DS = DynamicStallModels


airfoilpath = "./data/airfoils"
airfoil_interp_path = "./data/airfoils_interpolated"
yamlpath = "./data/5MW_PreComp_5seg"
precomppath = "./data/5MW_PreComp_5seg"
turbfile = "./data/TurbSim.dat"

# problem constants
rating = 5000.0 #Machine rating in kW (for cost model)
B = 3 #Number of blades
Rhub = 1.5 #meters, hub Radius
Rtip = 63.0 #meters, tip Radius
hubHt = 90.0 #meters, hub Height

precone = 2.5*pi/180 #radians
yaw = 0.0*pi/180 #radians
tilt = 0.0*pi/180 #radians
tsr0 = 7.55 #Initial tip speed ratio
pitch = 0.0 #parked pitch, radians

power_rating = 5.0e6 #Watts
max_thrust = 600e3 #Newtons
Vrated = 11.4 #m/s
Vinf = 10.0 #m/s
V_extreme = 70.0 #m/s
Vin = 3.0 #m/s
Vout = 25.0 #m/s
Vmean = 6.0 #m/s
Vtip = 80.0 #Tip speed in m/s (WISDEM example file)

shearExp = 0.2 #Shear exponent for wind profile
rho = 1.225 #kg/m^3, air density
mu = 1.837e-5 #kg/(m·s), dynamic viscosity of air
a = 335.0 #speed of sound, m/s
gravity = 9.81 #gravitational acceleration, m/s^2

azimuth0 = 0.0*pi/180 #initial azimuthal position, radians
azimuth = 90.0*pi/180 #extreme azimuth, radians

# tvec = collect(0:0.05:100.0) #time vector, seconds (for fatigue analysis) #Used 355 GB of RAM, too much for local testing.
@warn "Using shortened time vector for local testing."
tvec = collect(0:0.05:15.0)
ntime = length(tvec)
ntimecon = ntime - 200 #Number of constraints in time (dynamic tip deflection)

rotorR = Rtip*cos(precone)
Vcurve = collect(Vin:1.0:Vout) #wind speed vector for power/thrust curve
nwind = length(Vcurve)

### Read in the turbulent wind data
turb = readdlm(turbfile, skipstart=4)
n, m = size(turb)
tvec_turb = range(turb[1, 1], turb[n, 1], length=n) #Because the file doesn't save the time vector correctly. 
Ufit = Akima(tvec_turb, turb[:, 2])
Vfit = Akima(tvec_turb, turb[:, 5])
Wfit = Akima(tvec_turb, turb[:, 6])

env_data = (;Ufit, Vfit, Wfit)


L = Rtip - Rhub #Blade length
distro, materials, stations = uo.get_precomp_descriptions(precomppath, yamlpath)


m = 10.0 #Wholer exponent for fatigue
eps_ult = 0.01 #ultimate strain
nu = 1.3 #fatigue safety factor 
gamma_f = 1.35 #safety factor
gamma_m = 1.1 #safety factor
minf = 0.5674 #minimum thickness factor

nbearing = 2 #Number of bearings in the turbine

rvec = distro.rvec #meters, radial location of an analysis node
cvec = distro.cvec #meters, initial chord distribution
twistvec = distro.twistvec #Radians, initial twist distribution
le_loc = distro.le_loc #Location of reference axis as fraction of chord



### Create airfoil objects
nr = length(rvec)

airfoils = uo.get_interp_polars(rvec, airfoil_interp_path)


aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput(joinpath(airfoilpath, "Cylinder1.dat")) 
aftypes[2] = of.read_airfoilinput(joinpath(airfoilpath, "Cylinder2.dat")) 
aftypes[3] = of.read_airfoilinput(joinpath(airfoilpath, "DU40_A17.dat")) 
aftypes[4] = of.read_airfoilinput(joinpath(airfoilpath, "DU35_A17.dat")) 
aftypes[5] = of.read_airfoilinput(joinpath(airfoilpath, "DU30_A17.dat")) 
aftypes[6] = of.read_airfoilinput(joinpath(airfoilpath, "DU25_A17.dat")) 
aftypes[7] = of.read_airfoilinput(joinpath(airfoilpath, "DU21_A17.dat")) 
aftypes[8] = of.read_airfoilinput(joinpath(airfoilpath, "NACA64_A17.dat"))


raf = [2.8667, 5.6, 8.3333, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333]
afidx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]
af_names = ["Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat", "DU30_A17.dat", "DU25_A17.dat", "DU21_A17.dat", "NACA64_A17.dat"]

af_idx = of.integerfit(raf, afidx, rvec)

afs = aftypes[af_idx]


# dsairfoils = Vector{DS.Airfoil}(undef, nr) #Note: Apparently this can't display
dsairfoils = StructArray{DS.Airfoil}(undef, nr)
xcp = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    if isa(airfoils[i], CCBlade.Cylinder)
        polar_ = [-pi 0.0 airfoils[i].cd 0.0;
                 0.0 0.0 airfoils[i].cd 0.0;
                 pi 0.0 airfoils[i].cd 0.0]
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; dsmodel=DS.NoModel(), polar=polar_)
    else
        polar_ = hcat(airfoils[i].alpha, airfoils[i].cl, airfoils[i].cd, zeros(length(airfoils[i].alpha)))
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polar_)
    end
end



#Rotor object
rotor = Rotor(Rhub, Rtip, B, precone=precone, turbine=true)



### discretize the beam
pts = zeros(length(rvec)+1)
pts[1] = Rhub
for i = 1:length(rvec)-1
    pts[i+1] = (rvec[i] + rvec[i+1])/2
end
pts[end] = Rtip

points = [[pts[i], 0., 0.] for i in 1:length(pts)]
xp = [[rvec[i], 0., 0.0] for i in 1:length(rvec)]

nelem = length(points) - 1
start = 1:nelem
stop = 2:nelem+1

 




assembly = Assembly(points, start, stop; midpoints=xp)


# println("Functions scaling: ", string(Dates.now()))
### scaling factors
chord_scale = 1e1
twist_scale = 1e-1
thick_scale = 1e2
pitch_scale = 1e0
tsr_scale = 1e2


power_scale = 1e7
deflection_scale = 1e1
bending_scale = 2e7
buckling_scale = 1e4 #basically no change from 1e2 to 1e4...
strain_scale = 1e0
obj_scale = 1e2




num_elements = Int.(uo.num_elements5/2)
N_elements = Int(uo.N_elements5/2)

num_buckling = uo.num_buckling5
N_buckling = uo.N_buckling5

num_segments = uo.num_segments5 #Number of segments at each cross section
Nsegs = uo.Nsegs5

N_bulk_full = N_buckling + N_elements



fat_idxs = [1, 14, 16, 20, 22] #Locations to check for damage
Nfat = sum(Int, num_elements[fat_idxs])

cp_idxs = [1, 11, 18, 24, 32, nr] #Control point indices
nx_cp = length(cp_idxs) - 1

#14 is the end of the cylindrical section

twist_cp_idxs = [14, 18, 27, nr] #Control point indices for twist
nx_cp_twist = length(twist_cp_idxs) - 1

ntwist = nr - twist_cp_idxs[1] #Number of indices that we're controlling the twist at (minus 1). 
cyl_idxs = 1:twist_cp_idxs[1] #Indices of the cylindrical section
twist_idxs = twist_cp_idxs[1]:twist_cp_idxs[end] #Indices of the twist section



f_cp_idxs = [1, 14, 22, 31, 36]

nf_cp = length(f_cp_idxs)


function get_designvars(x, rvec, cvec, twistvec, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, nwind, chord_scale, twist_scale, thick_scale, tsr_scale, pitch_scale) 

    fit = (xx, yy) -> Akima(xx, yy, 1e-2)

    # x = x.*individual_scale #Scale the design variables by the individual scale. ->scaling turned off for derivative timing so we don't have to update the scaling parameters every time.

    ## Chords
    r_cp_chord = rvec[cp_idxs[1:end]] #Control point radii
    chord_idxs = 1:nx_cp
    x_chord = x[chord_idxs].*chord_scale
    chord_cp = vcat(cvec[cp_idxs[1]], x_chord) #Control point chords
    cfit = fit(r_cp_chord, chord_cp) 
    chords = cfit.(rvec) + zero(x[1:36])
 

    ## Twist
    start_idx = chord_idxs[end]
    twist_idxs = start_idx+1:start_idx+nx_cp_twist 
    r_cp_twist = rvec[twist_cp_idxs]
    x_twist = vcat(twistvec[twist_cp_idxs[1]], x[twist_idxs].*twist_scale)
    twistfit = fit(r_cp_twist, x_twist)
    cylinder_twists = twistvec[1:twist_cp_idxs[1]-1]
    blade_twists = twistfit.(rvec[twist_cp_idxs[1]:end])
    twists = vcat(cylinder_twists, blade_twists)


    for ti in x_twist 
        if ti >= pi/4
            println("Twist is too large.")
        end
    end
    

    ## segment scaling factors
    rsegs = rvec[f_cp_idxs]
    start_idx = twist_idxs[end]
    f1_idxs = start_idx+1:start_idx+nf_cp
    f1 = x[f1_idxs].*thick_scale #Scaling factors for the segments
    f1fit = fit(rsegs, f1)

    start_idx = f1_idxs[end]
    f2_idxs = start_idx+1:start_idx+nf_cp
    f2 = x[f2_idxs].*thick_scale
    f2fit = fit(rsegs, f2)

    start_idx = f2_idxs[end]
    f3_idxs = start_idx+1:start_idx+nf_cp
    f3 = x[f3_idxs].*thick_scale
    f3fit = fit(rsegs, f3)

    start_idx = f3_idxs[end]
    f4_idxs = start_idx+1:start_idx+nf_cp
    f4 = x[f4_idxs].*thick_scale
    f4fit = fit(rsegs, f4)

    start_idx = f4_idxs[end]
    f5_idxs = start_idx+1:start_idx+nf_cp
    f5 = x[f5_idxs].*thick_scale
    f5fit = fit(rsegs, f5)

    TF = typeof(x[1])
    fvec = zeros(TF, 5*length(rvec))
    for i in eachindex(rvec)
        idx = 5*(i-1) 
        fvec[idx+1] = f1fit(rvec[i])
        fvec[idx+2] = f2fit(rvec[i])
        fvec[idx+3] = f3fit(rvec[i])
        fvec[idx+4] = f4fit(rvec[i])
        fvec[idx+5] = f5fit(rvec[i])
    end

    ## pitches
    start_idx = f5_idxs[end]
    pitches_idxs = start_idx+1:start_idx+nwind
    pitches = x[pitches_idxs].*pitch_scale

    pitch0 = pitches[1]
    pitches = pitches .- pitch0 #Shift all the pitches by the first pitch. 
    twists = twists .+ pitch0 #Shift all the twists by the first pitch. 

    start_idx = pitches_idxs[end]
    tsr_idxs = start_idx+1
    tsr = x[tsr_idxs]*tsr_scale

    return chords, twists, fvec, pitches, tsr
end




aero = (; rotor, airfoils, dsairfoils, xcp, B, Rhub, Rtip, pitch, precone, yaw, tilt, azimuth0, azimuth, hubHt, shearExp, rho, mu, a, power_rating, max_thrust, Vrated, Vinf, Vin, Vout, V_extreme, Vmean, Vcurve, tsr0, env_data, nr, rotorR, nwind, rating, Vtip)

structural = (; xp, nelem, start, stop, points, tvec, gravity, nbearing)

composites = (; m, eps_ult, nu, gamma_f, gamma_m, minf, materials, stations)

idxs = (; fat_idxs, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, num_elements, num_buckling, N_elements, N_buckling, N_bulk_full, num_segments, Nsegs, Nfat, ntime, ntwist, cyl_idxs, twist_idxs, ntimecon)

scaling = (; chord_scale, thick_scale, twist_scale, power_scale, deflection_scale, bending_scale, obj_scale, buckling_scale, strain_scale, pitch_scale, tsr_scale)

objective = uo.ObjectiveFunction(distro, assembly, aero, structural, composites, idxs, scaling);
constraint = uo.ConstraintFunction(distro, assembly, aero, structural, composites, idxs, scaling);



function (obj::uo.ObjectiveFunction)(x; verbose::Bool=false, showfig::Bool=false) 

    TF = typeof(x[1])

    ### Read in design variables 
    chords, twists, fvec, pitches, tsr = get_designvars(
        x, 
        obj.distro.rvec, 
        obj.distro.cvec, 
        obj.distro.twistvec, 
        obj.idxs.cp_idxs, 
        obj.idxs.twist_cp_idxs, 
        obj.idxs.f_cp_idxs, 
        obj.idxs.nx_cp,
        obj.idxs.nx_cp_twist,
        obj.idxs.nf_cp,
        obj.aero.nwind,
        obj.scaling.chord_scale, 
        obj.scaling.twist_scale,
        obj.scaling.thick_scale,
        obj.scaling.tsr_scale,
        obj.scaling.pitch_scale)


    # println("Finished reading design variables.")
    #### Calculate the aerodynamic parts of the problem
    sections = CCBlade.Sections(obj.distro.rvec, chords, twists, obj.aero.airfoils)
    
    ### Power curve 
    _, _, AEP = uo.get_thrustpowercurve(
        obj.aero.rotor, 
        sections, 
        obj.aero.Vcurve, 
        tsr,
        pitches, 
        obj.aero.yaw, 
        obj.aero.tilt, 
        obj.aero.azimuth, 
        obj.aero.hubHt, 
        obj.aero.shearExp, 
        obj.aero.rho,
        obj.aero.Vmean,
        TF)

    ## Scale the layer thickness of the different sectors
    materials_i = GXBeamCS.Material{TF}.(obj.composites.materials) # Make appropriate materials vector. 
    segs_webs = uo.scale_segments(obj.composites.stations, fvec, materials_i, obj.idxs.num_segments)

    ### Get the cross sectional properties. 
    clt_sections = [uo.get_clt_sections_oop(
        obj.composites.stations[i].xaf, 
        obj.composites.stations[i].yaf, 
        chords[i], 
        twists[i], 
        obj.distro.le_loc[i], 
        obj.composites.stations[i].xbreak, 
        obj.composites.stations[i].webloc, 
        segs_webs[i].segments, 
        segs_webs[i].webs
    ) for i in 1:obj.aero.nr]

    mass_list = [GXBeamCS.mass_matrix(clt_sections[i]; reference=[chords[i]*obj.distro.le_loc[i], 0.0, twists[i]])[1] for i in eachindex(clt_sections)]

    ### Compute mass
    blade_mass = uo.compute_blade_mass(mass_list, obj.assembly.elements, TF)
    

    Q_rotor = uo.estimate_rotor_torque(obj.aero.rating, obj.aero.Rtip*2, tsr*obj.aero.Vinf) 

    tcc = uo.calc_turbine_cost(obj.aero.rating, obj.aero.Rtip*2, Q_rotor, obj.aero.B, blade_mass, obj.structural.nbearing, obj.aero.hubHt) #total capital cost

    bos = 2979e3 #Balance of Station cost (WISDEM NREL 5MW example file)
    opex = 144e3 #Annual O&M cost (WISDEM NREL 5MW example file)
    tr = 0.4 #Wisdem CSM default tax rate... 
    coe = uo.cost_of_energy(AEP/1e3, bos, tcc, opex; tax_rate=tr) #Convert the AEP from Wh to kWh

    if showfig
        @show AEP, blade_mass, coe
    end

    return obj.scaling.obj_scale*coe
end

function (obj::uo.ConstraintFunction)(g, x; verbose::Bool=false) 

    TF = typeof(x[1])

    ### Read in design variables 
    chords, twists, fvec, pitches, tsr = get_designvars(
        x, 
        obj.distro.rvec, 
        obj.distro.cvec, 
        obj.distro.twistvec, 
        obj.idxs.cp_idxs, 
        obj.idxs.twist_cp_idxs, 
        obj.idxs.f_cp_idxs, 
        obj.idxs.nx_cp,
        obj.idxs.nx_cp_twist,
        obj.idxs.nf_cp,
        obj.aero.nwind,
        obj.scaling.chord_scale, 
        obj.scaling.twist_scale,
        obj.scaling.thick_scale,
        obj.scaling.tsr_scale,
        obj.scaling.pitch_scale)

    ##### Geometric constraints
    ### Twist monotonicity constraint
    ntwist = obj.idxs.ntwist 
    twist_con_idxs = 1:ntwist #Indices of the twist constraints
    g[twist_con_idxs] = diff(twists[obj.idxs.twist_idxs]) #Checked for CA

    ### fvec constraint
    current_idx = twist_con_idxs[end]
    idx_fvec = current_idx+1:current_idx+length(fvec)
    g[idx_fvec] = fvec .- obj.composites.minf #Checked for CA


    #### Calculate the aerodynamic parts of the problem
    sections = CCBlade.Sections(obj.distro.rvec, chords, twists, obj.aero.airfoils)
    
    ### Power curve 
    thrustcurve, powercurve, AEP = uo.get_thrustpowercurve(
        obj.aero.rotor, 
        sections, 
        obj.aero.Vcurve, 
        tsr,
        pitches, 
        obj.aero.yaw, 
        obj.aero.tilt, 
        obj.aero.azimuth, 
        obj.aero.hubHt, 
        obj.aero.shearExp, 
        obj.aero.rho,
        obj.aero.Vmean,
        TF)

    ### Power constraints
    current_idx = idx_fvec[end] 
    power_idx1 = current_idx+1:current_idx+obj.aero.nwind 
    g[power_idx1] = -powercurve./obj.aero.power_rating #Greater than zero constraint 

    current_idx = power_idx1[end] 
    power_idx2 = current_idx+1:current_idx+obj.aero.nwind
    g[power_idx2] = powercurve./obj.aero.power_rating .- 1.0 #max power constraint. 


    ### Thrust constraints
    current_idx = power_idx2[end]
    thrust_idx = current_idx+1:current_idx+obj.aero.nwind
    g[thrust_idx] = thrustcurve./obj.aero.max_thrust .- 1.0 #max thrust constraint.


    ### pitches constraint
    current_idx_ = thrust_idx[end]
    idx_pitches = current_idx_+1:current_idx_+obj.aero.nwind-1
    g[idx_pitches] = -diff(pitches)
    

    ## Scale the layer thickness of the different sectors
    materials_i = GXBeamCS.Material{TF}.(obj.composites.materials) # Make appropriate materials vector. 
    segs_webs = uo.scale_segments(obj.composites.stations, fvec, materials_i, obj.idxs.num_segments) #TODO: Create new function to not scale fw1 and fw1

    
    ### Get the cross sectional properties. 
    closed_section = true

    clt_list = [GXBeamCS.CLT(GXBeamCS.get_beam_sections_oop(
        obj.composites.stations[i].xaf, 
        obj.composites.stations[i].yaf, 
        chords[i], 
        twists[i], 
        obj.distro.le_loc[i], 
        obj.composites.stations[i].xbreak, 
        obj.composites.stations[i].webloc, 
        segs_webs[i].segments, 
        segs_webs[i].webs), closed_section) for i in 1:obj.aero.nr]

    shear_center = true 
    

    compliance_list = [GXBeamCS.compliance_matrix(clt, shear_center)[1] for clt in clt_list]
    mass_list = [GXBeamCS.mass_matrix_clt(clt_list[i]; reference=[chords[i]*obj.distro.le_loc[i], 0.0, twists[i]])[1] for i in eachindex(clt_list)] 

    ## Calculate the aerodynamic loads at the extreme wind speed
    Omega = 0.0 #We should be in parked conditions... if the wind is at 70 m/s (156.6 mph), we should be parked.
    op = windturbine_op.(obj.aero.V_extreme, Omega, obj.aero.pitch, obj.distro.rvec, obj.aero.precone, obj.aero.yaw, obj.aero.tilt, obj.aero.azimuth, obj.aero.hubHt, obj.aero.shearExp, obj.aero.rho)
    out = CCBlade.solve.(Ref(obj.aero.rotor), sections, op)
    
    ### Interpolate the loads into the structural frame. #TODO: This comment probably needs some clarification. 
    fy, fz = uo.get_loads(out, chords, obj.aero.rho, TF)

    buckling, strain, static_tip_def = uo.extreme_loading_analysis_oop(
        clt_list,
        compliance_list,
        mass_list,
        fy,
        fz,
        obj.structural.points,
        obj.structural.xp,
        obj.aero.azimuth,
        obj.composites.eps_ult,
        obj.structural.gravity)

    ## buckling constraint
    current_idx = idx_pitches[end]
    buckling_idx = current_idx+1:current_idx+obj.idxs.N_buckling
    g[buckling_idx] = -buckling./obj.scaling.buckling_scale #Checked for CA

    ## Minimum strain constraint
    current_idx = buckling_idx[end]
    strain_idx1 = current_idx+1:current_idx+obj.idxs.N_elements
    g[strain_idx1] = -(strain./(obj.scaling.strain_scale*obj.composites.gamma_f*obj.composites.gamma_m) .+ 1) #Checked for CA

    ## Maximum strain constraint
    current_idx = strain_idx1[end]
    strain_idx2 = current_idx+1:current_idx+obj.idxs.N_elements
    g[strain_idx2] = strain./(obj.scaling.strain_scale*obj.composites.gamma_f*obj.composites.gamma_m) .- 1.0 #Checked for CA

    ## Static tip deflection constraint
    current_idx = strain_idx2[end]
    deflection_idx = current_idx+1
    g[deflection_idx] = 1.1 .* static_tip_def/obj.scaling.deflection_scale .+ 1


    Omega_rated = obj.aero.Vrated*tsr/obj.aero.rotorR 
    rated_pitch = pitches[9]
    dynamic_tip_deflections, damages = uo.fatigue_analysis_oop(
            clt_list,
            compliance_list,
            mass_list,
            chords,
            twists,
            Omega_rated,
            rated_pitch,
            obj.structural.points,
            obj.structural.xp,
            obj.structural.tvec,
            obj.aero.B,
            obj.distro.rvec,
            obj.aero.Rhub,
            obj.aero.Rtip,
            obj.aero.xcp,
            obj.aero.dsairfoils,
            obj.aero.hubHt,
            obj.aero.azimuth0,
            obj.aero.yaw,
            obj.aero.tilt,
            obj.aero.precone,
            obj.aero.env_data.Ufit,
            obj.aero.env_data.Vfit,
            obj.aero.env_data.Wfit,
            obj.aero.shearExp,
            obj.aero.rho,
            obj.aero.mu,
            obj.aero.a,
            obj.composites.eps_ult,
            obj.composites.m,
            obj.structural.gravity,
            obj.idxs.fat_idxs,
            obj.idxs.num_elements)

    current_idx = deflection_idx 
    tip_def_idx = current_idx+1:current_idx+obj.idxs.ntimecon
    g[tip_def_idx] = 1.1 .* dynamic_tip_deflections./obj.scaling.deflection_scale .+ 1  

    current_idx = tip_def_idx[end]
    damage_idxs = current_idx+1:current_idx+obj.idxs.Nfat
    g[damage_idxs] = log.(damages./obj.composites.nu)./200 #Dlife and scaling factor

end


### initial guess 
chords0 = cvec[cp_idxs[2:end]]./chord_scale
twist0 = twistvec[twist_cp_idxs[2:end]]./twist_scale
f_segs0 = ones(nf_cp*5)./thick_scale

pitches0 = zeros(nwind)./pitch_scale
pitches0[8:end] .= range(4.0*pi/180, stop=25.0*pi/180, length=nwind-7)./pitch_scale

tsr_naught = tsr0/tsr_scale


x0 = vcat(chords0, twist0, f_segs0, pitches0, tsr_naught)
x0 = x0

nx = length(x0)


### Bounds 
nc = length(chords0)
lb_chord = chords0.*(0.75/chord_scale) #This value is already scaled. 
ub_chord = (7.0/chord_scale).*ones(nc)

nt = length(twist0)
lb_twist = (0.0*pi/180).*ones(nt)./twist_scale 
ub_twist = (40.0*pi/180).*ones(nt)./twist_scale



lb_f_segs = minf*ones(5*nf_cp)./thick_scale
ub_f_segs = 3.0*ones(5*nf_cp)./thick_scale #5 segments at each cross section

lb_pitches = 0.0*pi/180*ones(nwind)./pitch_scale 
ub_pitches = 30.0*pi/180*ones(nwind)./pitch_scale

lb_tsr = 1.0/tsr_scale
ub_tsr = 15.0/tsr_scale

lx = vcat(lb_chord, lb_twist, lb_f_segs, lb_pitches, lb_tsr)
ux = vcat(ub_chord, ub_twist, ub_f_segs, ub_pitches, ub_tsr)

if length(lx) != length(ux) != nx
    @warn("Length of lower and upper bounds do not match the number of design variables.")
end

### Constraints
ng = ntwist + Nsegs + 3*nwind + (nwind-1) + N_buckling + 2*N_elements + 1 + ntimecon + Nfat #Number of constraints

lg = -Inf*ones(ng)
ug = zeros(ng)

constraint_names = ["Twist Monotonicity", "Fvec", "Power", "Thrust", "Pitches", "Buckling", "Strain", "Deflection", "Dynamic Tip Deflection", "Fatigue"]
constraint_nums = [ntwist, Nsegs, nwind, nwind, nwind-1, N_buckling, N_elements, 1, ntimecon, Nfat]

constraint_ends = cumsum(constraint_nums)





function make_cp_indices(start_idx::Int, end_idx::Int, n_cp::Int)
    # sanity check: enough room for n_cp distinct indices
    if n_cp > (end_idx - start_idx + 1)
        error("Cannot place $n_cp unique control points between $start_idx and $end_idx")
    end

    # create evenly spaced indices and round to Int
    idxs = round.(Int, range(start_idx, end_idx; length=n_cp))
    # ensure unique and sorted
    idxs = unique(sort(idxs))

    # if rounding caused duplicates, regenerate with a bit more resolution until we have enough
    while length(idxs) < n_cp
        candidate = round.(Int, range(start_idx, end_idx; length=length(idxs)+1))
        idxs = unique(sort(candidate))
    end

    return idxs
end


# nwind = parse(Int, ARGS[1]) #for use with slurm job arrays
nwind = 20
println("\nTesting with $nwind wind speeds...")

# Update wind speed curve and pitches
Vcurve = collect(Vin:((Vout-Vin)/(nwind-1)):Vout)
pitches0 = zeros(nwind)./pitch_scale
pitches0[8:end] .= range(4.0*pi/180, stop=25.0*pi/180, length=nwind-7)./pitch_scale 

# Redefine aero with new Vcurve and nwind
aero = (; rotor, airfoils, dsairfoils, xcp, B, Rhub, Rtip, pitch, precone, yaw, tilt, azimuth0, azimuth, hubHt, shearExp, rho, mu, a, power_rating, max_thrust, Vrated, Vinf, Vin, Vout, V_extreme, Vmean, Vcurve, tsr0, env_data, nr, rotorR, nwind, rating, Vtip)

# Initial guess for design variables
x0 = vcat(chords0, twist0, f_segs0, pitches0, tsr_naught)
x0 = x0 
nx = length(x0)
@show nx
flush(stdout)

# Redefine objective and constraint functions with new nwind
idxs = (; fat_idxs, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, num_elements, num_buckling, N_elements, N_buckling, N_bulk_full, num_segments, Nsegs, Nfat, ntime, ntwist, cyl_idxs, twist_idxs, ntimecon)
objective = uo.ObjectiveFunction(distro, assembly, aero, structural, composites, idxs, scaling)
constraint = uo.ConstraintFunction(distro, assembly, aero, structural, composites, idxs, scaling)

# Update constraint counts for new nwind
ng = ntwist + Nsegs + 3*nwind + (nwind-1) + N_buckling + 2*N_elements + 1 + ntimecon + Nfat

# Time function evaluation
println("Objective call...")
@time objective(deepcopy(x0))
@time objective(deepcopy(x0))
@time objective(deepcopy(x0))
@time objective(deepcopy(x0))

flush(stdout)

println("Constraint call...")
@time constraint(zeros(ng), deepcopy(x0))
# @time constraint(zeros(ng), deepcopy(x0))
# @time constraint(zeros(ng), deepcopy(x0))
# @time constraint(zeros(ng), deepcopy(x0))

flush(stdout)


function conwrap(x)
    TF = typeof(x[1])
    gg = zeros(TF, ng) #Initialize gg to be the same type as x for ReverseDiff compatibility.
    constraint(gg, x) 

    return ksmax(gg)
end

println("ReverseDiff....")
dgr = ReverseDiff.gradient(conwrap, deepcopy(x0))
@time ReverseDiff.gradient(conwrap, deepcopy(x0))
# @time ReverseDiff.gradient(conwrap, deepcopy(x0))
# @time ReverseDiff.gradient(conwrap, deepcopy(x0))
# @time ReverseDiff.gradient(conwrap, deepcopy(x0))


# println("Polyester....")
# dgp = zeros(ng)
# PolyesterForwardDiff.threaded_gradient!(conwrap, dgp, deepcopy(x0), ForwardDiff.Chunk(8))
# @time PolyesterForwardDiff.threaded_gradient!(conwrap, dgp, deepcopy(x0), ForwardDiff.Chunk(8))
# @time PolyesterForwardDiff.threaded_gradient!(conwrap, dgp, deepcopy(x0), ForwardDiff.Chunk(8))
# @time PolyesterForwardDiff.threaded_gradient!(conwrap, dgp, deepcopy(x0), ForwardDiff.Chunk(8))



nothing