#=
Parse the random start optimizations. 

Adam Cardoza
=# 
localpath = @__DIR__
cd(localpath)

using Pkg
Pkg.activate("./")

using Revise

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
rootname = "_COE_single_"


println("running ", filename, " at ", nowstr)

of = OpenFASTTools
uo = UnsteadyOpt
DS = DynamicStallModels


airfoilpath = "./data/Airfoils"
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

individual_scale = individual_scale .* Lscale 
# individual_scale_old = individual_scale .* Lscale 
# individual_scale = individual_scale .* Lscalenew
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
x0 = x0./individual_scale #Scale the initial guess by the individual scale.




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







constraint_names = ["Twist Monotonicity", "Fvec", "Power", "Thrust", "Pitches", "Buckling", "Strain", "Deflection", "Dynamic Tip Deflection", "Fatigue"]
constraint_nums = [ntwist, Nsegs, nwind, nwind, nwind-1, N_buckling, N_elements, 1, ntimecon, Nfat]

constraint_ends = cumsum(constraint_nums)



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
    fontfamily="times",
    color_palette=my_color_palette,
    grid=false,
    foreground_color_legend=nothing,
    tickfontsize=7,
    guidefontsize=7,
    legendfontsize=7,
)


# Parse all perturbation logs in current directory
if !@isdefined(parseflag)
    parseflag = true
end


if parseflag
    println("Parsing files.... ")
    all_results = uo.parse_all_perturbation_logs("./");
    parseflag = false
end


default(;
    # fontfamily="Palatino Roman",
    fontfamily="times",
    color_palette=my_color_palette,
    grid=false,
    foreground_color_legend=nothing,
    tickfontsize=7,
    guidefontsize=7,
    legendfontsize=7,
)

num_succeeded = count(res -> res.success, all_results)
num_converged = count(res -> res.converged, all_results)
num_failed = length(all_results) - num_succeeded
println("Finished optimizations: $num_succeeded")
println("Converged optimizations: $num_converged")
println("Failed optimizations: $num_failed")


failed_perturbations = [res.perturbation_num for res in all_results if !res.success]

chord_initalvals = []
twist_initialvals = []
pitch_initialvals = []
omega_initialvals = []

chord_vals = []
twist_vals = []
pitch_vals = []
omega_vals = [] 

div_group = Int[]
pernumbs = Int[]

cplt = plot(xaxis="Radius (m)", yaxis="Chord (m)", legend=true)
tplt = plot(xaxis="Radius (m)", yaxis="Twist (deg)", legend=false)
s1plt = plot(xaxis="Radius (m)", yaxis="Scaling Factor 1", legend=true)
s2plt = plot(xaxis="Radius (m)", yaxis="Scaling Factor 2", legend=true)
s3plt = plot(xaxis="Radius (m)", yaxis="Scaling Factor 3", legend=true)
s4plt = plot(xaxis="Radius (m)", yaxis="Scaling Factor 4", legend=true)
s5plt = plot(xaxis="Radius (m)", yaxis="Scaling Factor 5", legend=true)
pplt = plot(xaxis="Wind Speed (m/s)", yaxis="Pitch (deg)", legend=false)

labelbool = true
for (i, res) in enumerate(all_results)
    global labelbool
    if res.success
        chordsi, twistsi, fveci, pitchesi, tsri  = get_designvars(res.x0, rvec, cvec, twistvec, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, nwind, chord_scale, twist_scale, thick_scale, tsr_scale, pitch_scale, individual_scale)

        fmati = zeros(5, length(rvec))
        for i = 1:5
            fmati[i, :] = fveci[i:5:end]
        end

        chordso, twisto, fveco, pitcheso, tsro  = get_designvars(res.xopt, rvec, cvec, twistvec, cp_idxs, twist_cp_idxs, f_cp_idxs, nx_cp, nx_cp_twist, nf_cp, nwind, chord_scale, twist_scale, thick_scale, tsr_scale, pitch_scale, individual_scale)

        maxomega = 12.0*2*pi/60
        omega_initial = min.(Vcurve*tsri/rotorR, maxomega)
        omega_opt = min.(Vcurve*tsro/rotorR, maxomega)

        if labelbool
            plot!(cplt, rvec, chordsi, seriescolor=3, linealpha=0.5, lab="Initial Designs")
            plot!(tplt, rvec, twistsi.*(180/pi), seriescolor=3, linealpha=0.5, lab="Initial Designs")
            plot!(pplt, Vcurve, pitchesi.*(180/pi), seriescolor=3, linealpha=0.5, lab="Initial Designs")
        else
            plot!(cplt, rvec, chordsi, seriescolor=3, linealpha=0.5, lab=false)
            plot!(tplt, rvec, twistsi.*(180/pi), seriescolor=3, linealpha=0.5, lab=false)
            plot!(pplt, Vcurve, pitchesi.*(180/pi), seriescolor=3, linealpha=0.5, lab=false)
        end

        push!(chord_initalvals, chordsi)
        push!(twist_initialvals, twistsi)
        push!(pitch_initialvals, pitchesi)
        push!(omega_initialvals, omega_initial)

        push!(chord_vals, chordso)
        push!(twist_vals, twisto)
        push!(pitch_vals, pitcheso)
        push!(omega_vals, omega_opt)
        push!(pernumbs, res.perturbation_num)

        fmato = zeros(5, length(rvec))
        for i = 1:5
            fmato[i, :] = fveco[i:5:end]
        end

        if maximum(chordso) > 6
            push!(div_group, res.perturbation_num)
        end

        if labelbool
            plot!(cplt, rvec, chordso, seriescolor=1, lab="Final Designs", lw=1.5)
            plot!(tplt, rvec, twisto.*(180/pi), seriescolor=1, lab="Final Designs", lw=0.5)
            plot!(pplt, Vcurve, pitcheso.*(180/pi), seriescolor=1, lab="Final Designs", lw=0.5)
        else
            plot!(cplt, rvec, chordso, seriescolor=1, lab=false, lw=1.5)
            plot!(tplt, rvec, twisto.*(180/pi), seriescolor=1, lab=false, lw=1.5)
            plot!(pplt, Vcurve, pitcheso.*(180/pi), seriescolor=1, lab=false, lw=1.5)
        end
        labelbool=false
    end
end

using Statistics

# Convert chord_vals to a matrix where each column is a chord distribution
chord_matrixi = hcat(chord_initalvals...)
twist_matrixi = hcat(twist_initialvals...)
pitch_matrixi = hcat(pitch_initialvals...)
omega_matrixi = hcat(omega_initialvals...)

chord_matrix = hcat(chord_vals...)
twist_matrix = hcat(twist_vals...)
pitch_matrix = hcat(pitch_vals...)
omega_matrix = hcat(omega_vals...)

# using JLD2
# @save "randomstart_variable.jld2" rvec Vcurve chord_matrixi twist_matrixi pitch_matrixi omega_matrixi chord_matrix twist_matrix pitch_matrix omega_matrix

# Calculate mean and std for each radius location (row-wise)
chord_mean = mean(chord_matrix, dims=2)
chord_median = median(chord_matrix, dims=2)
chord_std = std(chord_matrix, dims=2)

twist_mean = mean(twist_matrix, dims=2)
twist_median = median(twist_matrix, dims=2)
twist_std = std(twist_matrix, dims=2)

pitch_mean = mean(pitch_matrix, dims=2)
pitch_median = median(pitch_matrix, dims=2)
pitch_std = std(pitch_matrix, dims=2)




display(cplt)
display(tplt)
display(pplt)


function in2figsize(width_in, height_in)
    # px_per_inch = 96  # PGFplots default
    px_per_inch = 72 # GR default
    width_px = width_in * px_per_inch
    height_px = height_in * px_per_inch
    return (round(Int, width_px), round(Int, height_px))
end

# textwidth = 7.2 #WES
textwidth = 6.5 #AIAA

ratio = (8, 4)
width = 0.95*textwidth
height = width*ratio[2]/ratio[1]
figsize = in2figsize(width, height)



layout = @layout [ a{0.6w} [grid(2,1)] ]

randplt = plot(cplt, tplt, pplt, layout=layout, size=figsize, dpi=300)
display(randplt)
# savefig(randplt, "randomopt_all_"*nowstr*".pdf")



nothing