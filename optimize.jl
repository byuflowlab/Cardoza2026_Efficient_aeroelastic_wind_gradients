#=
Optimize a simplified version of the NREL 5MW wind turbine blade. 
objective: COE
by varying: chord, twist, thickness scaling factors, pitch schedule, tip speed ratio
constraints: power, thrust, deflection, strains, buckling, fatigue damage #Todo: 

Adam Cardoza 8/20/25
=# 
localpath = @__DIR__
cd(localpath)

using GXBeamCS, GXBeam, CCBlade, OpenFASTTools, DynamicStallModels
using FLOWMath, DelimitedFiles, LinearAlgebra, Dates
using StaticArrays, StructArrays, SparseArrays
using ForwardDiff, FiniteDiff, DiffResults, PolyesterForwardDiff
using SNOW, Snopt
using Plots, LaTeXStrings
using UnsteadyOpt

df = Dates.DateFormat("yymmdd_HH.MM.SS")
now = Dates.now()
nowstr = Dates.format(now, df)
filename = splitpath(@__FILE__)[end]
rootname = "_COE_opt_"


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

tvec = collect(0:0.05:100.0) #time vector, seconds (for fatigue analysis)
# tvec = collect(0:0.025:100.0) 
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


L = Rtip - Rhub #blade length
distro, materials, stations = uo.get_precomp_descriptions(precomppath, yamlpath)


m = 10.0 #Wholer exponent
eps_ult = 0.01 #Ultimate strain for all materials (assumed)
nu = 1.3 #fatigue failure safety factor
nbearing = 2 #Number of bearings in the turbine

rvec = distro.rvec #meters, radial location of an analysis node
cvec = distro.cvec #meters, initial chord distribution
twistvec = distro.twistvec #Radians, initial twist distribution
le_loc = distro.le_loc #Location of reference axis as fraction of chord


### Create airfoil objects
nr = length(rvec)

airfoils = uo.get_interp_polars(rvec, airfoil_interp_path) #Interpolated polars (no dynamic coefficients)


aftypes = Array{of.AirfoilInput}(undef, 8) #OpenFAST airfoils (have dynamic coefficients), but not interpolated to station
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



### scaling factors for optimization
chord_scale = 1e1
twist_scale = 1e-1
thick_scale = 1e2
pitch_scale = 1e0
tsr_scale = 1e2


power_scale = 1e7
thrust_scale = 600e3 #600 kN max thrust constraint. 
deflection_scale = 1e1
bending_scale = 2e7
buckling_scale = 1e4 #basically no change from 1e2 to 1e4...
strain_scale = 1e0
obj_scale = 1e2

individual_scale = [10.0, 10.0, 1.0, 1.0, 1.0, 10.0, 1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 100.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 0.1, 0.1, 1.0, 10.0, 1.0, 1.0, 10.0, 10.0, 100.0, 1000.0, 100.0, 100.0, 100.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 1.0] #Scaling to get the gradient to order zero. 

Lscale = [0.1, 0.1, 1.0, 10.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 1.0, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 1.0, 0.1, 0.01, 10.0, 10.0, 1.0, 0.1, 0.01, 0.1, 1.0, 0.01, 0.1, 0.01, 0.001, 0.1, 1.0, 1.0, 1.0, 100.0, 0.1, 1000.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0, 100.0, 10.0, 100.0, 100.0, 0.1] #Lagrangian based-scaling

Lscalenew = [10.0, 10.0, 10.0, 0.1, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 10.0, 10.0, 10.0, 0.1, 1.0, 10.0, 100.0, 10.0, 1.0, 10.0, 1.0, 10.0, 100.0, 0.1, 1.0, 10.0, 10.0, 100.0, 10.0, 1.0, 1.0, 1.0, 100.0, 100.0, 0.1, 1.0, 0.1, 1.0, 0.01, 10.0, 0.01, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 10.0] 

# individual_scale = individual_scale .* Lscale 
# individual_scale_old = individual_scale .* Lscale 
individual_scale = individual_scale .* Lscalenew

@warn "Using Lagrangian scaled inputs"


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

f_cp_idxs = [1, 14, 22, 31, 36] #control point indices for thickness scaling factors
nf_cp = length(f_cp_idxs)


function get_designvars(x, rvec, cvec, twistvec, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, nwind, chord_scale, twist_scale, thick_scale, tsr_scale, pitch_scale, individual_scale) 
    fit = (xx, yy) -> Akima(xx, yy, 1e-4)

    x = x.*individual_scale #Scale the design variables by the individual scale.

    ## Chords
    r_cp_chord = rvec[cp_idxs[1:end]] #Control point radii
    chord_idxs = 1:nx_cp
    x_chord = x[chord_idxs].*chord_scale
    chord_cp = vcat(cvec[cp_idxs[1]], x_chord) #Control point chords
    cfit = fit(r_cp_chord, chord_cp) 
    chords = cfit.(rvec)

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




aero = (; rotor, airfoils, dsairfoils, xcp, B, Rhub, Rtip, pitch, precone, yaw, tilt, azimuth0, azimuth, hubHt, shearExp, rho, mu, a, power_rating, Vrated, Vinf, Vin, Vout, V_extreme, Vmean, Vcurve, tsr0, env_data, nr, rotorR, nwind, rating, Vtip)

structural = (; xp, nelem, start, stop, points, tvec, gravity, nbearing)

composites = (; m, eps_ult, nu, materials, stations)

idxs = (; fat_idxs, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, num_elements, num_buckling, N_elements, N_buckling, N_bulk_full, num_segments, Nsegs, Nfat, ntime, ntwist, cyl_idxs, twist_idxs, ntimecon)

scaling = (; chord_scale, thick_scale, twist_scale, power_scale, thrust_scale, deflection_scale, bending_scale, obj_scale, buckling_scale, strain_scale, pitch_scale, tsr_scale, individual_scale)

objective = uo.ObjectiveFunction(distro, assembly, aero, structural, composites, idxs, scaling);
constraint = uo.ConstraintFunction(distro, assembly, aero, structural, composites, idxs, scaling);



function (obj::uo.ObjectiveFunction)(x; verbose::Bool=false) 

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
        obj.scaling.pitch_scale,
        obj.scaling.individual_scale)


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
        TF) #searched for branching. 


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
    
    tcc = uo.calc_turbine_cost(obj.aero.rating, obj.aero.Rtip*2, Q_rotor, obj.aero.B, blade_mass, obj.structural.nbearing, obj.aero.hubHt)

    bos = 2979e3 #Balance of Station cost (WISDEM NREL 5MW example file)
    opex = 144e3 #Annual O&M cost (WISDEM NREL 5MW example file)
    tr = 0.4 #Wisdem CSM default tax rate
    coe = uo.cost_of_energy(AEP/1e3, bos, tcc, opex; tax_rate=tr) #Convert the AEP from Wh to kWh

    if verbose
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
        obj.scaling.pitch_scale,
        obj.scaling.individual_scale)


    ##### Geometric constraints
    ### Twist monotonicity constraint
    ntwist = obj.idxs.ntwist 
    twist_con_idxs = 1:ntwist #Indices of the twist constraints
    g[twist_con_idxs] = diff(twists[obj.idxs.twist_idxs]) 


    ### fvec constraint
    current_idx = twist_con_idxs[end]
    idx_fvec = current_idx+1:current_idx+length(fvec)
    g[idx_fvec] = fvec 


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
    power_idx = current_idx+1:current_idx+obj.aero.nwind 
    g[power_idx] = powercurve./obj.scaling.power_scale 

    ### Thrust constraints
    current_idx = power_idx[end]
    thrust_idx = current_idx+1:current_idx+obj.aero.nwind
    g[thrust_idx] = thrustcurve./obj.scaling.thrust_scale


    ### pitches constraint
    current_idx_ = thrust_idx[end]
    idx_pitches = current_idx_+1:current_idx_+nwind-1
    g[idx_pitches] = diff(pitches)
    

    ## Scale the layer thickness of the different sectors
    materials_i = GXBeamCS.Material{TF}.(obj.composites.materials) # Make appropriate materials vector. 
    segs_webs = uo.scale_segments(obj.composites.stations, fvec, materials_i, obj.idxs.num_segments) 

    
    ### Get the cross sectional properties. 
    clt_list = [uo.get_clt(
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

    shear_center = true 

    compliance_list = [GXBeamCS.compliance_matrix(clt, shear_center)[1] for clt in clt_list]
    mass_list = [GXBeamCS.mass_matrix_clt(clt_list[i]; reference=[chords[i]*obj.distro.le_loc[i], 0.0, twists[i]])[1] for i in eachindex(clt_list)] 



    ## Calculate the aerodynamic loads at the extreme wind speed
    Omega = 0.0 #parked conditions
    op = windturbine_op.(obj.aero.V_extreme, Omega, obj.aero.pitch, obj.distro.rvec, obj.aero.precone, obj.aero.yaw, obj.aero.tilt, obj.aero.azimuth, obj.aero.hubHt, obj.aero.shearExp, obj.aero.rho)
    out = CCBlade.solve.(Ref(obj.aero.rotor), sections, op)
    
    ### Dimensionalize the loads and rotate into GXBeam reference frame
    fy, fz = uo.get_loads(out, chords, obj.aero.rho, TF)

    current_idx = idx_pitches[end]
    uo.extreme_loading_analysis!(
        g,
        clt_list,
        compliance_list,
        mass_list,
        fy,
        fz,
        obj.structural.points,
        obj.structural.xp,
        obj.aero.azimuth,
        obj.composites.eps_ult,
        obj.structural.gravity,
        obj.idxs.num_elements,
        obj.idxs.num_buckling,
        obj.idxs.N_buckling,
        obj.idxs.N_elements,
        current_idx,
        obj.scaling.buckling_scale,
        obj.scaling.strain_scale,
        obj.scaling.deflection_scale)

    
    current_idx = current_idx + obj.idxs.N_elements + obj.idxs.N_buckling + 1 #Update the current index to the end of the extreme loading analysis.
    Omega_rated = obj.aero.Vrated*tsr/obj.aero.rotorR 
    rated_pitch = pitches[9]
    uo.fatigue_analysis!(
            g,
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
            obj.composites.nu,
            obj.structural.gravity,
            obj.idxs.fat_idxs,
            obj.idxs.num_elements,
            obj.idxs.N_bulk_full,
            current_idx,
            obj.idxs.ntimecon,
            obj.scaling.deflection_scale)
    

end

### initial guess 
chords0 = cvec[cp_idxs[2:end]]./chord_scale
twist0 = twistvec[twist_cp_idxs[2:end]]./twist_scale
f_segs0 = ones(nf_cp*5)./thick_scale

pitches0 = zeros(nwind)./pitch_scale
pitches0[8:end] .= range(4.0*pi/180, stop=25.0*pi/180, length=nwind-7)./pitch_scale 

tsr_naught = tsr0/tsr_scale


x0 = vcat(chords0, twist0, f_segs0, pitches0, tsr_naught)
x0 = x0./individual_scale #Scale the initial guess by the individual scale. -> Auto includes the new scaling 

# x0 = [0.4028239358431749, 0.4675537947166109, 0.4010315053735982, 0.02899559634967371, 0.1586684149921422, 0.16448182148645088, 0.6610865363107865, 0.01457692416100335, 0.009066056025266573, 0.010252321258681528, 0.00951991951516333, 0.00947690246977771, 0.010172193690986481, 0.010415441473851268, 0.10351677043895527, 0.10706151335039121, 0.10057912915908272, 0.010267724625617767, 0.010271759406807999, 0.10942056110269935, 0.01006434098877854, 0.0960372277723652, 0.09514015116970374, 0.0009558701923253718, 0.009678358260949843, 0.1066419185056164, 0.10801140177499911, 0.09374282121096847, 0.10148177894656758, 0.010291509798503674, 0.10889338351082843, 0.010925879643683316, 0.01074552195837471, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.279137596645605e-5, 0.08625786993578798, 0.12190877911997504, 0.1463207664646107, 0.15666498517184385, 0.18970625219078685, 0.19728098283592002, 0.24841191564895287, 0.02840684031478236, 0.028131748159335244, 0.028875893241993206, 0.030830941425237356, 0.0033467136958755445, 0.003965520489397802, 0.00039814475828889534, 0.0004070550004539656, 0.8078374272442163] #Using old scaling #random start that has a negative under the square root. 

dv_lengths = [length(chords0), length(twist0), length(f_segs0), length(pitches0), 1]
dv_names = ["Chords", "Twists", "Segment Thickness Scaling Factors", "Pitches", "Tip Speed Ratio"]
dv_ends = cumsum(dv_lengths)

nx = length(x0)



### Bounds 
nc = length(chords0)
lb_chord = chords0.*(0.75/chord_scale) #This value is already scaled. 
ub_chord = (7.0/chord_scale).*ones(nc)

nt = length(twist0)
lb_twist = (0.0*pi/180).*ones(nt)./twist_scale 
ub_twist = (40.0*pi/180).*ones(nt)./twist_scale

minf = 0.5674 

lb_f_segs = minf*ones(5*nf_cp)./thick_scale
ub_f_segs = 3.0*ones(5*nf_cp)./thick_scale #5 segments at each cross section

lb_pitches = 0.0*pi/180*ones(nwind)./pitch_scale
ub_pitches = 30.0*pi/180*ones(nwind)./pitch_scale

lb_tsr = 1.0/tsr_scale
ub_tsr = 15.0/tsr_scale

lx = vcat(lb_chord, lb_twist, lb_f_segs, lb_pitches, lb_tsr)./individual_scale
ux = vcat(ub_chord, ub_twist, ub_f_segs, ub_pitches, ub_tsr)./individual_scale

if length(lx) != length(ux) != nx
    @warn("Length of lower and upper bounds do not match the number of design variables.")
end

### Constraints
ng = ntwist + Nsegs + 2*nwind + (nwind-1) + N_buckling + N_elements + 1 + ntimecon + Nfat #Number of constraints

## Geometry constraints
lg_twist_monotonicity = -Inf.*ones(ntwist)
ug_twist_monotonicity = zeros(ntwist)

lg_fvec = minf.*ones(Nsegs)
ug_fvec = 3.0.*ones(Nsegs)

## power case constraints
lg_power = zeros(nwind)
ug_power = ones(nwind).*(5e6/power_scale) 

lg_thrust = zeros(nwind)
ug_thrust =  ones(nwind)

lg_pitches = zeros(nwind-1)
ug_pitches = ones(nwind-1).*Inf 

## Extreme case constraints
lg_buckling = zeros(N_buckling)
ug_buckling = Inf.*ones(N_buckling)     

gamma_f = 1.35 #safety factor
gamma_m = 1.1
lg_strain = -ones(N_elements)./(gamma_f*gamma_m)
ug_strain = ones(N_elements)./(gamma_f*gamma_m) 

lg_deflection = -5.0191*0.9/deflection_scale
ug_deflection = abs(lg_deflection) 

### Unsteady case constraints
lg_dyn_tip_deflection = (-5.0191*0.9/deflection_scale).*ones(ntimecon)
ug_dyn_tip_deflection = abs.(lg_dyn_tip_deflection)

lg_fatigue = -Inf.*ones(Nfat)
ug_fatigue = zeros(Nfat) # log(1) = 0 => constrains lifetime damage to less than 1.0. 


lg = vcat(lg_twist_monotonicity, lg_fvec, lg_power, lg_thrust, lg_pitches, lg_buckling, lg_strain, lg_deflection, lg_dyn_tip_deflection, lg_fatigue)
ug = vcat(ug_twist_monotonicity, ug_fvec, ug_power, ug_thrust, ug_pitches, ug_buckling, ug_strain, ug_deflection, ug_dyn_tip_deflection, ug_fatigue)

length(lg) != length(ug) != ng ? @warn("Length of lower and upper bounds do not match the number of constraints.") : nothing


# @time f0 = objcon(g0, deepcopy(x0); showfig=false) 
# @time f0 = objcon(g0, deepcopy(x0); showfig=false) 
# @time f0 = objcon(g0, deepcopy(x0); showfig=false) 






constraint_names = ["Twist Monotonicity", "Fvec", "Power", "Thrust", "Pitches", "Buckling", "Strain", "Deflection", "Dynamic Tip Deflection", "Fatigue"]
constraint_nums = [ntwist, Nsegs, nwind, nwind, nwind-1, N_buckling, N_elements, 1, ntimecon, Nfat]

constraint_ends = cumsum(constraint_nums)


# f0 = objective(deepcopy(x0); verbose=true) #Initial objective function value

# g0 = zeros(ng) #Initialize g0 to zero.
# constraint(g0, deepcopy(x0); verbose=false) #Precompile the constraint function (and so we can check initial condition)

# near_tol = 9e-6  # You can adjust this threshold
# near_violations = uo.get_near_violations(g0, lg, ug, near_tol, constraint_names, constraint_ends)
# println("Active constraints:")
# display(near_violations)
# println("")








using Statistics


derivative_option = UserDeriv()
rc = readdlm("./data/sparsity_pattern.csv", ',', Int)
sparsity_pattern = sparse(rc[:,1], rc[:,2], ones(length(rc[:,1])), ng, nx) #Sparsity pattern for the jacobian. 
chunksize = 8
deriv = uo.SparseForwardDeriv(objective, constraint, nx, ng, sparsity_pattern, chunksize, deepcopy(x0))

snopt_options = Dict("Major iterations limit" => 600,
# "Minor iterations limit" => 10,
"Time Limit" => 45*3600, # 45 hours
"Derivative option" => 1, #Derivatives are known
"Verify level" => 0, #Check the derivatives
"Major feasibility tolerance" => 9e-6,
"Major optimality tolerance" => 5e-5, 
"Print file" => "snopt-print"*filename*"_"*rootname*"_"*nowstr*".out",
"Summary file" => "snopt-summary"*filename*"_"*rootname*"_"*nowstr*".out",
"New basis file" => 17,
"Save frequency" => 100)

# options = Options(solver=IPOPT(), derivatives=derivative_option)
options = Options(solver=SNOPT(options=snopt_options), derivatives=derivative_option)

run_optimization = false
if run_optimization
    println("Starting $(splitpath(@__FILE__)[end]) optimization...")
    flush(stdout)
    
    xopt, fopt, info, out = minimize(deriv, deepcopy(x0), ng, lx, ux, lg, ug, options)


    println("Optimization complete.")
    println("x0: ", x0)
    println("f0: ", f0)
    println("xopt: ", xopt)
    println("fopt: ", fopt)
    println("info: ", info)

else
    # xopt = [0.4651205911237532, 0.5775436523147818, 0.32828553967715285, 0.021466763443071853, 0.09350195616919597, 0.1693408336024107, 0.5944998820790519, 0.0, 0.0056739813319321665, 0.00644937180954058, 0.0102716537689427, 0.014257629236724668, 0.00652508538502962, 0.005688139964188984, 0.061212235399601744, 0.08651341318646076, 0.11476646633862035, 0.0062237605427874765, 0.005674, 0.09954084903561385, 0.01285377157048517, 0.2100037596074625, 0.1200357886177187, 0.0005674002571620446, 0.005674027680370433, 0.05675027614266712, 0.09810956473727252, 0.07844561326368246, 0.05693542036678655, 0.006313164268127356, 0.09965019274353461, 0.014141119246335862, 0.006003330873326401, 0.0, 0.0, 0.0, 0.0, 0.0, 2.8019301190203405e-21, 0.0, -1.1828532072194644e-22, 0.01630541121142068, 0.06225645448334745, 0.11265368819973096, 0.1518221984405436, 0.1856560160919358, 0.2162040233554313, 0.2444401499969176, 0.027081679313868203, 0.029564605252868338, 0.03194504025161199, 0.03423415116165135, 0.0036425098131472374, 0.0038536768092994464, 0.00040655203914800945, 0.00043746524732403576, 0.899504474082129] #Lagrangian scale optimization #Using old scaling. 

    # xopt = [0.4593766274842783, 0.6079654678708539, 0.34745173898124254, 0.021799630523988967, 0.09761840083350544, 0.17700838617687184, 0.6220324969867851, 0.0, 0.00567662261831312, 0.006266891978110481, 0.009434461707465834, 0.014255938199946234, 0.005716303023396108, 0.00569391626946602, 0.06038927058980086, 0.08373663657410874, 0.1160930896593463, 0.007180872556696759, 0.0056739999999999985, 0.10048804765362664, 0.011038527646940154, 0.2023088037383343, 0.11902092790509292, 0.0005673999999999999, 0.005674312255854106, 0.05675171562003144, 0.09339560663434233, 0.07933552218007471, 0.056959445283055035, 0.0062675424716054084, 0.09727946593871861, 0.01379198359286497, 0.006672718938919495, 0.0, 1.6163557037957021e-15, 1.6158339315001935e-16, 1.6158339315001935e-16, 1.6163557037957021e-15, 1.6138095227562557e-17, 1.6139867220488213e-14, 1.6147836106455982e-17, 0.01932420182331069, 0.06446167931331649, 0.11397058311798172, 0.15268914699257702, 0.18620685700207967, 0.21650528091119287, 0.24453263616352733, 0.027073216217443272, 0.02954208686708262, 0.031909858919764174, 0.03418764021016623, 0.003636912387662408, 0.0038471893229299367, 0.00040512987795814873, 0.00042499975807188127, 0.8712594162729468] #Near optimal solution for unsteady fix. #Using old scaling

    xopt = [0.004639483227113929, 0.005933254440389638, 0.03452180428033821, 2.180160186292858, 0.09763364541378818, 0.17612391962821178, 0.6212305085173769, 0.0, 0.000569435248393383, 0.0006317173554264929, 9.956409193967199e-5, 0.0001403018583133886, 6.6228733365260324e-6, 0.05687657912535405, 0.0061031954727681975, 0.0008534157723675139, 0.00011261541057711974, 6.689763998865057e-5, 0.005673999999999992, 0.000999318953196919, 0.011639541205320958, 0.002008334527756839, 1.2076337736613505e-5, 0.05673999999999996, 0.05674290380441979, 0.005675089451049552, 0.0009636514977382998, 7.756421564581263e-6, 0.0005693794324622878, 0.006296347477088405, 0.000982099669351496, 0.001378041061803274, 6.372943266452945e-7, 0.0, 0.0, 6.777695796028599e-21, 6.777695634469885e-20, 0.0, 0.0, -1.5812354499253308e-20, 0.0, 0.0018708794879829526, 0.06396286113628272, 0.11349653857454063, 0.1522232491762067, 0.18574154768241122, 0.21603678453013744, 0.02440603156894168, 0.27025950710422864, 0.2949395638633526, 0.318599952677618, 0.03413731256694071, 0.36354708293167126, 0.03841995966241092, 0.040459656906239226, 0.04244525195041182, 0.00871000231741218] #Optimal solution for unsteady fix. #Using new scaling.

    # xopt = [0.4028239358431749, 0.4675537947166109, 0.4010315053735982, 0.02899559634967371, 0.1586684149921422, 0.16448182148645088, 0.6610865363107865, 0.01457692416100335, 0.009066056025266573, 0.010252321258681528, 0.00951991951516333, 0.00947690246977771, 0.010172193690986481, 0.010415441473851268, 0.10351677043895527, 0.10706151335039121, 0.10057912915908272, 0.010267724625617767, 0.010271759406807999, 0.10942056110269935, 0.01006434098877854, 0.0960372277723652, 0.09514015116970374, 0.0009558701923253718, 0.009678358260949843, 0.1066419185056164, 0.10801140177499911, 0.09374282121096847, 0.10148177894656758, 0.010291509798503674, 0.10889338351082843, 0.010925879643683316, 0.01074552195837471, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.279137596645605e-5, 0.08625786993578798, 0.12190877911997504, 0.1463207664646107, 0.15666498517184385, 0.18970625219078685, 0.19728098283592002, 0.24841191564895287, 0.02840684031478236, 0.028131748159335244, 0.028875893241993206, 0.030830941425237356, 0.0033467136958755445, 0.003965520489397802, 0.00039814475828889534, 0.0004070550004539656, 0.8078374272442163] #Using old scaling #random start that has a negative under the square root. 
end





run_analysis = false
if run_analysis

    println("Running optimization result analysis...")

    gopt = zeros(ng)
    println("   Running objective function...")
    ffinal = objective(deepcopy(xopt); verbose=true)
    println("   Running constraint function...")
    constraint(gopt, deepcopy(xopt); verbose=false)

    gopt_violations = uo.check_constraints(gopt, lg, ug)

    near_tol = 5e-5  # You can adjust this threshold
    near_constraints_opt = uo.get_near_violations(gopt, lg, ug, near_tol, constraint_names, constraint_ends) 

    @show near_constraints_opt

    y = zeros(ng)
    x = deepcopy(xopt)
    Jopt = zeros(ng, nx) #Jacobian matrix


    println("   Running sparse derivative...")
    fderiv = uo.SparseForwardDeriv(objective, constraint, nx, ng, sparsity_pattern, chunksize, deepcopy(x0))
    gopt3 = zeros(ng)
    df_opt3 = zeros(nx)
    dg_opt3 = zeros(ng, nx)
    f_opt3 = fderiv(gopt3, df_opt3, dg_opt3, deepcopy(xopt))

    println("   Running Finite difference derivatives...")
    finderiv = uo.FiniteDeriv(objective, constraint, nx, ng, sparsity_pattern)
    gopt4 = zeros(ng)
    df_opt4 = zeros(nx)
    dg_opt4 = zeros(ng, nx)
    f_opt4 = finderiv(gopt4, df_opt4, dg_opt4, deepcopy(xopt))


    ### Get the scaling factors based on the Jacobian
    # Jscale = uo.scale_jacobian(dg_opt3; smax=15)

    ### Get the scaling factors based on the gradient of the Lagrangian. 
    # fmul = reshape(readdlm("./data/coe_mm_fmul.csv", ','), :)
    # xmul = reshape(readdlm("./data/coe_mm_xmul.csv", ','), :) 

    # dLdx = @. df_opt3 + xmul

    # for i = 1:ng
    #     dLdx .+= @. fmul[i] * dg_opt3[i, :]
    # end

    # Lscale = uo.scale_gradient(dLdx; smax=15)
end



############## Extract optimization results and plot
plot_results = true #Runs the plots code (Rewrites some variables above)
show_fig = false
save_fig = false
base_name = "_unsteady_extrachord_rotated_"
if @isdefined(xopt) && plot_results
    using Plots.Measures

    my_color_palette=[
        RGB(0.0, 46.0 / 255.0, 93.0 / 255.0), #BYU Blue
        RGB(155.0 / 255.0, 0.0, 0.0), #"BYU" Red
        RGB(128.0 / 255.0, 128.0 / 255.0, 128.0 / 255.0), #Middle Gray
        RGB(162.0 / 255.0, 227.0 / 255.0, 162.0 / 255.0), #Light Green
        RGB(243.0 / 255.0, 209.0 / 255.0, 243.0 / 255.0), #Pink
        RGB(205.0 / 255.0, 179.0 / 255.0, 0.0), #Yellow
        RGB(161.0 / 255.0, 161.0 / 255.0, 226.0 / 255.0), #Purple
    ]
    default(;
        fontfamily="Palatino Roman",
        color_palette=my_color_palette,
        grid=false,
        foreground_color_legend=nothing,
    )
    fit = (xx, yy) -> Akima(xx, yy, 1e-4)

    using Colors
    theme_colors = palette(:auto)
    filename = splitpath(@__FILE__)[end]


    chords0, twists0, fvec0, pitches_initial, tsr_initial  = get_designvars(x0, rvec, cvec, twistvec, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, nwind, chord_scale, twist_scale, thick_scale, tsr_scale, pitch_scale, individual_scale)
    chords_opt, twists_opt, fvec_opt, pitches_opt, tsr_opt  = get_designvars(xopt, rvec, cvec, twistvec, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, nwind, chord_scale, twist_scale, thick_scale, tsr_scale, pitch_scale, individual_scale)

    rsmooth = range(rvec[1], rvec[end], length=100)
    rsmooth_webs = range(rvec[7], rvec[end], length=100)

    xopt_plot = xopt .* individual_scale

    r_cp = rvec[cp_idxs[1:end]] 
    chord_idxs = 1:nx_cp
    x_chord = xopt_plot[chord_idxs].*scaling.chord_scale
    c_cp = vcat(cvec[cp_idxs[1]], x_chord) 
    chords0fit = fit(rvec, cvec)
    chords_optfit = fit(rvec, chords_opt)
    chords0_smooth = chords0fit.(rsmooth)
    chords_opt_smooth = chords_optfit.(rsmooth)


    start_idx = chord_idxs[end]
    twist_idxs = start_idx+1:start_idx+nx_cp_twist 
    r_cp_twist = rvec[twist_cp_idxs[2:end]]
    x_twist = xopt_plot[twist_idxs].*scaling.twist_scale
    twists0fit = fit(rvec, twistvec)
    twists_optfit = fit(rvec, twists_opt)
    twists0_smooth = twists0fit.(rsmooth)
    twists_opt_smooth = twists_optfit.(rsmooth)

    start_idx = twist_idxs[end]
    seg_idxs = start_idx+1:start_idx+nf_cp
    f1_cp_opt = xopt_plot[seg_idxs].*thick_scale
    r_f_cp = rvec[f_cp_idxs]

    start_idx = seg_idxs[end]
    seg_idxs = start_idx+1:start_idx+nf_cp
    f2_cp_opt = xopt_plot[seg_idxs].*thick_scale

    start_idx = seg_idxs[end]
    seg_idxs = start_idx+1:start_idx+nf_cp
    f3_cp_opt = xopt_plot[seg_idxs].*thick_scale

    start_idx = seg_idxs[end]
    seg_idxs = start_idx+1:start_idx+nf_cp
    f4_cp_opt = xopt_plot[seg_idxs].*thick_scale

    start_idx = seg_idxs[end]
    seg_idxs = start_idx+1:start_idx+nf_cp
    f5_cp_opt = xopt_plot[seg_idxs].*thick_scale

    fmat0 = zeros(5, length(rvec))
    fmat_opt = zeros(5, length(rvec))
    for i = 1:5
        fmat0[i, :] = fvec0[i:5:end]
        fmat_opt[i, :] = fvec_opt[i:5:end]
    end



    seg_thickness0 = uo.get_seg_thickness(objective, chords0, twists0, fvec0) 
    seg_thicknessopt = uo.get_seg_thickness(objective, chords_opt, twists_opt, fvec_opt) #Doesn't actually need the objective, just the parameters

    maxomega = 12.0*2*pi/60
    omega_initial = min.(Vcurve*tsr_initial/rotorR, maxomega)
    omega_opt = min.(Vcurve*tsr_opt/rotorR, maxomega)




    c1plt = plot(rsmooth, chords0_smooth, xlabel="Radius (m)", ylabel="Chord (m)", lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, ylims=(0, 8))
    plot!(c1plt, rsmooth, chords_opt_smooth, lab="Optimized", lw=2)
    scatter!(c1plt, r_cp, c_cp, lab=false, seriescolor=2)
    # show_fig ? display(c1plt) : nothing
    # save_fig ? savefig(c1plt, filename*base_name*"chord_"*nowstr*".png") : nothing

    # c2plt = plot(objective, cvec, aspect_ratio=:equal, legend=false, xaxis=false, yaxis=false, grid=false, fillcolor=1)
    c2plt = plot(objective, chords_opt, aspect_ratio=:equal, legend=false, xaxis=false, yaxis=false, grid=false, fillcolor=2, xlims=xlims(c1plt), fillalpha=0.8)

    l = @layout [a{0.875h}; b{0.125h}]
    cplt = plot(c1plt, c2plt, layout=l, size=(800, 750))
    # display(cplt)

    twistplt = plot(rsmooth, twists0_smooth, xlabel="Radius (m)", ylabel="Twist (rad)", lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, size=(800, 600), legend=false)
    plot!(twistplt, rsmooth, twists_opt_smooth, lab="Optimized", lw=2)
    scatter!(twistplt, r_cp_twist, x_twist, lab=false, seriescolor=2)
    # show_fig ? display(twistplt) : nothing
    # save_fig ? savefig(twistplt, filename*base_name*"twist_"*nowstr*".png") : nothing

    pitchplt = plot(Vcurve, pitches_initial, xlabel="Wind Speed (m/s)", ylabel="Pitch (rad)", lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, size=(800, 600), legend=false)
    plot!(pitchplt, Vcurve, pitches_opt, lab="Optimized", lw=2)
    # show_fig ? display(pitchplt) : nothing
    # save_fig ? savefig(pitchplt, filename*base_name*"pitch_"*nowstr*".png") : nothing
    ax2 = twinx()
    plot!(ax2, Vcurve, omega_initial, lab="Original", ylabel="Rotational Speed (rad/s)", lw=2, legend=false, seriescolor=1, linestyle=:dash)
    plot!(ax2, Vcurve, omega_opt, lab="Optimized", lw=2, legend=false, seriescolor=2, linestyle=:dash)

    


    scalingplt = plot(rvec, fmat_opt[1, :], xlabel="Radius (m)", ylabel="Optimized Scaling Factor", lab="Seg 1", lw=2, seriescolor=1, grid=false, foreground_color_legend=nothing, background_color_legend=nothing)
    for i in 2:5
        if i < 4
            coloridx = i
        else
            coloridx = i + 2
        end
        plot!(scalingplt, rvec, fmat_opt[i, :], lab="Seg $i", lw=2, seriescolor=coloridx)
    end
    
    scatter!(scalingplt, r_f_cp, f1_cp_opt, lab=false, seriescolor=1)
    scatter!(scalingplt, r_f_cp, f2_cp_opt, lab=false, seriescolor=2)
    scatter!(scalingplt, r_f_cp, f3_cp_opt, lab=false, seriescolor=3)
    scatter!(scalingplt, r_f_cp, f4_cp_opt, lab=false, seriescolor=6)
    scatter!(scalingplt, r_f_cp, f5_cp_opt, lab=false, seriescolor=7)

    # show_fig ? display(scalingplt) : nothing
    # save_fig ? savefig(scalingplt, filename*base_name*"scaling_"*nowstr*".png") : nothing

    # thickplt = plot(rvec, seg_thicknessopt[:,1:5], xlabel="Radius (m)", ylabel="Segment Thickness (m)", lab=[L"S_1" L"S_2" L"S_3" L"S_4" L"S_5"], lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing)
    # display(thickplt)

    thick1plt = plot(rvec, seg_thickness0[:,1].*1e3, ylabel=L"S_1",  lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, bordercolor=my_color_palette[3], axiscolor=my_color_palette[3], fg_color_text=my_color_palette[3], bottom_margin=-3.25mm, legend=false)
    plot!(thick1plt, rvec, seg_thicknessopt[:,1].*1e3, lab="Optimized", lw=2)
    # display(thick1plt)

    thick2plt = plot(rvec, seg_thickness0[:,2].*1e3, ylabel=L"S_2",  lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[4], axiscolor=my_color_palette[4], fg_color_text=my_color_palette[4], top_margin=-2.5mm, bottom_margin=-2.5mm)
    plot!(thick2plt, rvec, seg_thicknessopt[:,2].*1e3, lab="Optimized", lw=2)
    # display(thick2plt)

    thick3plt = plot(rvec, seg_thickness0[:,3].*1e3, ylabel=L"S_3",  lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[5], axiscolor=my_color_palette[5], fg_color_text=my_color_palette[5], top_margin=-2.5mm, bottom_margin=-2.5mm)
    plot!(thick3plt, rvec, seg_thicknessopt[:,3].*1e3, lab="Optimized", lw=2)
    annotate!(thick3plt, -14.0, 30, text("Segment Thickness (mm)", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(thick3plt)

    thick4plt = plot(rvec, seg_thickness0[:,4].*1e3, ylabel=L"S_4",  lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[6], axiscolor=my_color_palette[6], fg_color_text=my_color_palette[6], top_margin=-2.5mm, bottom_margin=-2.5mm)
    plot!(thick4plt, rvec, seg_thicknessopt[:,4].*1e3, lab="Optimized", lw=2)
    # display(thick4plt)

    thick5plt = plot(rvec, seg_thickness0[:,5].*1e3, xlabel="Radius (m)", ylabel=L"S_5",  lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[7], axiscolor=my_color_palette[7], fg_color_text=my_color_palette[7], top_margin=-2.5mm, ylims=(0, 80))
    plot!(thick5plt, rvec, seg_thicknessopt[:,5].*1e3, lab="Optimized", lw=2)
    # display(thick5plt)
    # savefig(thick5plt, filename*base_name*"thickness_seg5_"*nowstr*".png")

    thickplt = plot(thick1plt, thick2plt, thick3plt, thick4plt, thick5plt, layout=(5, 1), size=(800, 600), left_margin=10mm, bordercolor=:black, axiscolor=:black, fg_color_text=:black)
    # display(thickplt)

    pltidx = 22
    crossplt = plot(objective, chords0, twists_opt, fvec_opt, pltidx, legend=false, xaxis=false, yaxis=false, grid=false)
    annotate!(crossplt, -0.05, -0.65, text(L"S_1", :left, 10, :black))
    annotate!(crossplt, 0.55, -0.65, text(L"S_2", :left, 10, :black))
    annotate!(crossplt, 1.55, -0.65, text(L"S_3", :left, 10, :black))
    annotate!(crossplt, 2.9, -0.65, text(L"S_4", :left, 10, :black))
    annotate!(crossplt, 3.9, -0.65, text(L"S_5", :left, 10, :black))
    display(crossplt)

    l2 = @layout [a{0.875h}; b{0.125h}]
    # layupplt = plot(scalingplt, crossplt, layout=l2, size=(800, 750))
    layupplt = plot(thickplt, crossplt, layout=l2, size=(800, 750))

    crossfig = plot(crossplt, aspect_ratio=:equal, size=(800,300))
    # savefig(crossfig, filename*base_name*"cross_section_"*nowstr*".png")
    ### scaling factor on top right
    # lf = @layout [a{0.6h} b{0.6h}; c{0.4h} d{0.4h}]
    # optplt = plot(cplt, layupplt, twistplt, pitchplt, layout=lf)
    # display(optplt)
    # savefig(optplt, filename*base_name*"optimization_results_"*nowstr*".png")



    #### Scaling factors on bottom
    angleplt = plot(twistplt, pitchplt, layout=(2, 1), size=(800, 600), legend=false)

    l3 = @layout [a{0.875h}; b{0.125h}]
    toprightplt = plot(angleplt, crossplt, layout=l3, size=(800, 750))


    cylidxs = 1:14
    cyl1plt = plot(rvec[cylidxs], seg_thickness0[cylidxs,1].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, bottom_margin=-2.0mm, legend=false, yticks=myticks(3, seg_thickness0[cylidxs,1].*1e3, seg_thicknessopt[cylidxs,1].*1e3))
    plot!(cyl1plt, rvec[cylidxs], seg_thicknessopt[cylidxs,1].*1e3, lab="Optimized", lw=2)
    annotate!(cyl1plt, -0.1, 32, text(L"S_1", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(cyl1plt)

    cyl2plt = plot(rvec[cylidxs], seg_thickness0[cylidxs,2].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[4], axiscolor=my_color_palette[4], fg_color_text=my_color_palette[4], top_margin=-2.5mm, bottom_margin=-2.5mm, yticks=myticks(3, seg_thickness0[cylidxs,2].*1e3, seg_thicknessopt[cylidxs,2].*1e3))
    plot!(cyl2plt, rvec[cylidxs], seg_thicknessopt[cylidxs,2].*1e3, lab="Optimized", lw=2)
    annotate!(cyl2plt, -0.1, 33, text(L"S_2", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(cyl2plt)

    cyl3plt = plot(rvec[cylidxs], seg_thickness0[cylidxs,3].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[5], axiscolor=my_color_palette[5], fg_color_text=my_color_palette[5], top_margin=-2.5mm, bottom_margin=-2.5mm, yticks=myticks(3, seg_thickness0[cylidxs,3].*1e3, seg_thicknessopt[cylidxs,3].*1e3))
    plot!(cyl3plt, rvec[cylidxs], seg_thicknessopt[cylidxs,3].*1e3, lab="Optimized", lw=2)
    annotate!(cyl3plt, -0.1, 42, text(L"S_3", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    annotate!(cyl3plt, -1.6, 30, text("Cylinder", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    annotate!(cyl3plt, -0.9, 30, text("Segment Thickness (mm)", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(cyl3plt)

    cyl4plt = plot(rvec[cylidxs], seg_thickness0[cylidxs,4].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[6], axiscolor=my_color_palette[6], fg_color_text=my_color_palette[6], top_margin=-2.5mm, bottom_margin=-2.5mm, yticks=myticks(3, seg_thickness0[cylidxs,4].*1e3, seg_thicknessopt[cylidxs,4].*1e3))
    plot!(cyl4plt, rvec[cylidxs], seg_thicknessopt[cylidxs,4].*1e3, lab="Optimized", lw=2)
    annotate!(cyl4plt, -0.1, 64, text(L"S_4", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(cyl4plt)

    cyl5plt = plot(rvec[cylidxs], seg_thickness0[cylidxs,5].*1e3, xlabel="Radius (m)", lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, bordercolor=my_color_palette[7], axiscolor=my_color_palette[7], fg_color_text=my_color_palette[7], top_margin=-2.5mm, yticks=myticks(3, seg_thickness0[cylidxs,5].*1e3, seg_thicknessopt[cylidxs,5].*1e3))
    plot!(cyl5plt, rvec[cylidxs], seg_thicknessopt[cylidxs,5].*1e3, lab="Optimized", lw=2)
    annotate!(cyl5plt, -0.1, 58, text(L"S_5", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(cyl5plt)
    # savefig(cyl5plt, filename*base_name*"thickness_seg5_"*nowstr*".png")

    cylplt = plot(cyl1plt, cyl2plt, cyl3plt, cyl4plt, cyl5plt, layout=(5, 1), size=(800, 600), left_margin=16mm, bordercolor=:black, axiscolor=:black, fg_color_text=:black)
    # display(cylplt)


    bladeidxs = 15:nr
    blade1plt = plot(rvec[bladeidxs], seg_thickness0[bladeidxs,1].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, bottom_margin=-2.0mm, legend=false, yticks=myticks(3, seg_thickness0[bladeidxs,1].*1e3, seg_thicknessopt[bladeidxs,1].*1e3), ylims=(2.85, 8.5))
    plot!(blade1plt, rvec[bladeidxs], seg_thicknessopt[bladeidxs,1].*1e3, lab="Optimized", lw=2)
    annotate!(blade1plt, 5, 6, text(L"S_1", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(blade1plt)

    blade2plt = plot(rvec[bladeidxs], seg_thickness0[bladeidxs,2].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, top_margin=-2.5mm, bottom_margin=-2.5mm, yticks=myticks(3, seg_thickness0[bladeidxs,2].*1e3, seg_thicknessopt[bladeidxs,2].*1e3), ylims=(5, 30.670875651062577))
    plot!(blade2plt, rvec[bladeidxs], seg_thicknessopt[bladeidxs,2].*1e3, lab="Optimized", lw=2)
    annotate!(blade2plt, 5, 17, text(L"S_2", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(blade2plt)

    blade3plt = plot(rvec[bladeidxs], seg_thickness0[bladeidxs,3].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, top_margin=-2.5mm, bottom_margin=-2.5mm, yticks=myticks(3, seg_thickness0[bladeidxs,3].*1e3, seg_thicknessopt[bladeidxs,3].*1e3))
    plot!(blade3plt, rvec[bladeidxs], seg_thicknessopt[bladeidxs,3].*1e3, lab="Optimized", lw=2)
    annotate!(blade3plt, 5, 37, text(L"S_3", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    annotate!(blade3plt, 1, 30, text("Segment Thickness (mm)", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    annotate!(blade3plt, -2.5, 30, text("Airfoil", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(blade3plt)

    blade4plt = plot(rvec[bladeidxs], seg_thickness0[bladeidxs,4].*1e3, lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, fg_color_text=my_color_palette[6], top_margin=-2.5mm, bottom_margin=-2.5mm, yticks=myticks(3, seg_thickness0[bladeidxs,4].*1e3, seg_thicknessopt[bladeidxs,4].*1e3))
    plot!(blade4plt, rvec[bladeidxs], seg_thicknessopt[bladeidxs,4].*1e3, lab="Optimized", lw=2)
    annotate!(blade4plt, 5, 53, text(L"S_4", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(blade4plt)

    blade5plt = plot(rvec[bladeidxs], seg_thickness0[bladeidxs,5].*1e3, xlabel="Radius (m)", lab="Original", lw=2, grid=false, foreground_color_legend=nothing, background_color_legend=nothing, legend=false, fg_color_text=my_color_palette[7], top_margin=-2.5mm, yticks=myticks(3, seg_thickness0[bladeidxs,5].*1e3, seg_thicknessopt[bladeidxs,5].*1e3), ylims=(2.8, 9.5))
    plot!(blade5plt, rvec[bladeidxs], seg_thicknessopt[bladeidxs,5].*1e3, lab="Optimized", lw=2)
    annotate!(blade5plt, 5, 6, text(L"S_5", :black, font(11, "Palatino Roman"), :center, rotation = 90))
    # display(blade5plt)
    # savefig(blade5plt, filename*base_name*"thickness_seg5_"*nowstr*".png")

    bladeplt = plot(blade1plt, blade2plt, blade3plt, blade4plt, blade5plt, layout=(5, 1), size=(800, 600), left_margin=10mm, bordercolor=:black, axiscolor=:black, fg_color_text=:black)
    # display(bladeplt)


    lf = @layout [a{0.6h} b{0.6h}; c{0.4h} d{0.4h}]
    optplt = plot(cplt, toprightplt, cylplt, bladeplt, layout=lf)
    display(optplt)
    # savefig(optplt, filename*base_name*"optimization_results_focusthicknesses_"*nowstr*".png")
    # savefig(optplt, filename*base_name*"optimization_results_focusthicknesses_"*nowstr*".pdf")
    nothing
end

nothing