



using OpenFASTTools, DelimitedFiles, GXBeam, WATT, LinearAlgebra, DynamicStallModels
using StaticArrays, StructArrays, Statistics
using Plots, Plots.Measures, LaTeXStrings
using Dates, FLOWMath



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
guidefontsize=7,
legendfontsize=7,
tickfontsize=7,
)

df = Dates.DateFormat("yymmdd_HH.MM.SS")
now = Dates.now()
nowstr = Dates.format(now, df)



DS = DynamicStallModels
of = OpenFASTTools

localpath = @__DIR__
cd(localpath)

using Dates
df = Dates.DateFormat("yymmdd_HH.MM.SS")
nowtime = Dates.now()
nowstr = Dates.format(nowtime, df)

fileending = ".png"

turbfile = "./data/TurbSim.dat"



### Read in OpenFAST files
ofpath = "./data/" 
inputfile = of.read_inputfile("sn5_input.fst", ofpath)
inflowwind = of.read_inflowwind("sn5_inflowwind.dat", ofpath)
# addriver = of.read_addriver("sn5_ADdriver.dvr", ofpath)
# adfile = of.read_adfile("sn5_ADfile.dat", ofpath)
adblade = of.read_adblade("sn5_ADblade.dat", ofpath)
edfile = of.read_edfile("sn5_edfile.dat", ofpath)
bdfile = of.read_bdfile("sn5_bdfile.dat", ofpath)
bdblade = of.read_bdblade("sn5_bdblade.dat", ofpath)

B = edfile["NumBl"] #Number of Blades

rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]
rfrac = bdblade["rfrac"]
chordvec = adblade["BlChord"]
twistvec = adblade["BlTwist"]

### Simulation corrections? 
rvec[1] += 0.001 
rtip = 63.0

hubht = inflowwind["RefHt"]
pitch = edfile["BlPitch(1)"]
precone = edfile["PreCone(1)"]*pi/180 
yaw = edfile["NacYaw"]*(pi/180) 
tilt = edfile["ShftTilt"]*(pi/180)

rho = inputfile["AirDens"]
mu = inputfile["KinVisc"]*rho
vinf = inflowwind["HWindSpeed"]
a = inputfile["SpdSound"]
omega = edfile["RotSpeed"]*2*pi/60 #Convert to rads/s
shearexp = inflowwind["PLExp"]
shearexp = 0.1
tsr = omega*rtip/vinf


turb = readdlm(turbfile, skipstart=4)
n, m = size(turb)
tvec_turb = range(turb[1, 1], turb[n, 1], length=n) #Because the file doesn't get the time vector correctly. 
Ufit = Akima(tvec_turb, turb[:, 2])
Vfit = Akima(tvec_turb, turb[:, 5])
Wfit = Akima(tvec_turb, turb[:, 6])

ufun(t) = SVector(Ufit(t), Vfit(t), Wfit(t))
omegafun(t) = SVector(0.0, 0.0, 0.0)
udotfun(t) = SVector(0.0, 0.0, 0.0) #Note: Could probably use the gradient function. 
omegadotfun(t) = SVector(0.0, 0.0, 0.0)
Vinf(t) = sqrt(Ufit(t)^2 + Vfit(t)^2)
RS(t) = omega
Vinfdot(t) = 0.0
RSdot(t) = 0.0
env = WATT.SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)




n = length(adblade["BlSpn"])
ne = Int(bdblade["station_total"])

if !@isdefined(resetflag)
    resetflag = false
end

if resetflag
    readflag = true
    defflag = true
    runflag = true

    resetflag = false
end

if !@isdefined(readflag)
    readflag = true
end

if readflag
    println("Reading OpenFAST files...")
    fullouts = readdlm("./data/sn5_input.out", skipstart=6)

    names = fullouts[1,:]

    data = Float64.(fullouts[3:end,:])

    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]
    


    nt = length(tvec)


    global fxmat = zeros(nt, n)
    global fymat = zeros(nt, n)
    global Mmat = zeros(nt, n)

    for i = 1:n
        if i<10
            number = "00$i"
        elseif i<100
            number = "0$i"
        else
            number = "$i"
        end
        namex = "AB1N"*number*"Fx"
        namey = "AB1N"*number*"Fy"
        nameux = "AB1N"*number*"Vx" 
        nameuy = "AB1N"*number*"Vy"
        namem = "AB1N"*number*"Mm"

        fxmat[:,i] = outs[namex]
        fymat[:,i] = outs[namey]
        Mmat[:,i] = outs[namem]
    end

    readflag = false
end

azimuth = outs["Azimuth"]
tipdx = outs["B1TipTDxr"]
tipdy = outs["B1TipTDyr"]
tipdz = outs["B1TipTDzr"]


assembly = of.make_assembly(edfile, bdfile, bdblade)

### Prep the ASD rotor and operating conditions 
aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput("./data/Airfoils/Cylinder1.dat") 
aftypes[2] = of.read_airfoilinput("./data/Airfoils/Cylinder2.dat") 
aftypes[3] = of.read_airfoilinput("./data/Airfoils/DU40_A17.dat") 
aftypes[4] = of.read_airfoilinput("./data/Airfoils/DU35_A17.dat") 
aftypes[5] = of.read_airfoilinput("./data/Airfoils/DU30_A17.dat") 
aftypes[6] = of.read_airfoilinput("./data/Airfoils/DU25_A17.dat") 
aftypes[7] = of.read_airfoilinput("./data/Airfoils/DU21_A17.dat") 
aftypes[8] = of.read_airfoilinput("./data/Airfoils/NACA64_A17.dat") 

# indices correspond to which airfoil is used at which station
af_idx = Int.(adblade["BlAFID"])


# create airfoil array
afs = aftypes[af_idx]

n = length(rvec)
airfoils = StructArray{DS.Airfoil}(undef, n)
xcp = Vector{Float64}(undef, n)
for i = 1:n
    airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
end

rR = rvec./rtip
blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; rhub=rhub, rtip=rtip, precone)

turbine = true
rotor_r = WATT.Rotor(Int(B), hubht, turbine; tilt, yaw)



if !@isdefined(defflag)
    defflag = true
end


if defflag
    println("initialize...")
    aerostates, gxhistory, mesh = WATT.initialize_sim(blade, assembly, tvec; verbose=true)
    defflag = false
end

if !@isdefined(runflag)
    runflag = true
end


if runflag 
    println("Running simulation...")
    println(Dates.now())
    WATT.run_sim!(rotor_r, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false)
    println(Dates.now())

    runflag = false 
end

# @time WATT.run_sim!(rotor_r, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false)
# @time WATT.run_sim!(rotor_r, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false)
# @time WATT.run_sim!(rotor_r, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false)
# @time WATT.run_sim!(rotor_r, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false)
# @time WATT.run_sim!(rotor_r, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false)


#Tip deflections
ntr = length(tvec)
tipdef_x = [gxhistory[i].points[end].u[1] for i in eachindex(tvec)]
tipdef_y = [gxhistory[i].points[end].u[2] for i in eachindex(tvec)]
tipdef_z = [gxhistory[i].points[end].u[3] for i in eachindex(tvec)]

tiptheta_x = zeros(ntr)
tiptheta_y = zeros(ntr)
tiptheta_z = zeros(ntr)

tiptheta_xof = zeros(nt)
tiptheta_yof = zeros(nt)
tiptheta_zof = zeros(nt)

for i = 1:ntr
    theta = WATT.WMPtoangle(gxhistory[i].points[end].theta)
    tiptheta_x[i] = theta[1]
    tiptheta_y[i] = theta[2]
    tiptheta_z[i] = theta[3]
end

for i = 1:nt #nt
    thetawmp = SVector(outs["B1TipRDxr"][i], outs["B1TipRDyr"][i], outs["B1TipRDzr"][i])
    theta = WATT.WMPtoangle(thetawmp)
    tiptheta_xof[i] = theta[1]
    tiptheta_yof[i] = theta[2]
    tiptheta_zof[i] = theta[3]
end

nr = length(rvec)
nt = length(tvec)
defx_gx = zeros(ntr, nr)
defy_gx = zeros(ntr, nr)
defz_gx = zeros(ntr, nr)

for i in eachindex(tvec)
    for j in eachindex(rvec)
        defx_gx[i,j] = gxhistory[i].points[j].u[1]
        defy_gx[i,j] = gxhistory[i].points[j].u[2]
        defz_gx[i,j] = gxhistory[i].points[j].u[3]
    end
end




using Plots, LaTeXStrings

function calc_error(a, b)
    return @. abs(a - b)
end

function calc_relative_error(a, b)
    return @. 100*(a - b)/b
end

function calc_normalized_error(a, b)
    normfactor = maximum(abs.(b))
    # @show normfactor
    return @. 100*(a - b)/(sign(b)*normfactor)
end


# idxs = 1:50:length(tvec)
idxs = 1:length(tvec)
## dt = 0.01
# testidxs = 501:length(tvec) #5-100s
testidxs = 1001:length(tvec) #10-100s
# testidxs = 5001:length(tvec) #50-100s

## dt = 0.05
# testidxs = 101:length(tvec) #5-100s
# testidxs = 201:length(tvec) #10-100ss
# testidxs = 1001:length(tvec) #50-100s

## dt = 0.001
# testidxs = 5001:length(tvec) #5-100s
# testidxs = 10001:length(tvec) #10-100s
# testidxs = 50001:length(tvec) #50-100s

fxerr_rel = calc_relative_error(aerostates.Fx[idxs,end][testidxs], fxmat[testidxs,end])
fyerr_rel = calc_relative_error(aerostates.Fy[idxs,end][testidxs], -fymat[testidxs,end])
fxerr = calc_error(aerostates.Fx[idxs,end][testidxs], fxmat[testidxs,end])
fyerr = calc_error(aerostates.Fy[idxs,end][testidxs], -fymat[testidxs,end])
fxerr_norm = calc_normalized_error(aerostates.Fx[idxs,end][testidxs], fxmat[testidxs,end])
fyerr_norm = calc_normalized_error(aerostates.Fy[idxs,end][testidxs], -fymat[testidxs,end])

tdxerr_rel = calc_relative_error( -tipdef_z[idxs][testidxs], tipdx[testidxs])
tdyerr_rel = calc_relative_error( tipdef_y[idxs][testidxs], tipdy[testidxs])
tdzerr_rel = calc_relative_error( tipdef_x[idxs][testidxs], tipdz[testidxs])
tdxerr = calc_error( -tipdef_z[idxs][testidxs], tipdx[testidxs])
tdyerr = calc_error( tipdef_y[idxs][testidxs], tipdy[testidxs])
tdzerr = calc_error( tipdef_x[idxs][testidxs], tipdz[testidxs])
tdxerr_norm = calc_normalized_error( -tipdef_z[idxs][testidxs], tipdx[testidxs])
tdyerr_norm = calc_normalized_error( tipdef_y[idxs][testidxs], tipdy[testidxs])
tdzerr_norm = calc_normalized_error( tipdef_x[idxs][testidxs], tipdz[testidxs])

println("Time averaged relative Tip Loading error (%): ")
@show maximum(abs.(fxerr_rel)), mean(abs.(fxerr_rel)), std(abs.(fxerr_rel))
@show maximum(abs.(fyerr_rel)), mean(abs.(fyerr_rel)), std(abs.(fyerr_rel))
println("Time averaged relative Tip deflection error (%): ")
@show maximum(abs.(tdxerr_rel)), mean(abs.(tdxerr_rel)), std(abs.(tdxerr_rel))
@show maximum(abs.(tdyerr_rel)), mean(abs.(tdyerr_rel)), std(abs.(tdyerr_rel))
@show maximum(abs.(tdzerr_rel)), mean(abs.(tdzerr_rel)), std(abs.(tdzerr_rel))
println("")


println("Time averaged absolute Tip Loading error: ")
@show maximum(abs.(fxerr)), mean(abs.(fxerr)), std(abs.(fxerr))
@show maximum(abs.(fyerr)), mean(abs.(fyerr)), std(abs.(fyerr))

println("Time averaged absolute Tip deflection error: ")
@show maximum(abs.(tdxerr)), mean(abs.(tdxerr)), std(abs.(tdxerr))
@show maximum(abs.(tdyerr)), mean(abs.(tdyerr)), std(abs.(tdyerr))
@show maximum(abs.(tdzerr)), mean(abs.(tdzerr)), std(abs.(tdzerr))

println("")
println("Time averaged normalized Tip Loading error: ")
@show maximum(abs.(fxerr_norm)), mean(abs.(fxerr_norm)), std(abs.(fxerr_norm))
@show maximum(abs.(fyerr_norm)), mean(abs.(fyerr_norm)), std(abs.(fyerr_norm))
println("Time averaged normalized Tip deflection error: ")
@show maximum(abs.(tdxerr_norm)), mean(abs.(tdxerr_norm)), std(abs.(tdxerr_norm))
@show maximum(abs.(tdyerr_norm)), mean(abs.(tdyerr_norm)), std(abs.(tdyerr_norm))
@show maximum(abs.(tdzerr_norm)), mean(abs.(tdzerr_norm)), std(abs.(tdzerr_norm))



idxs = 1:length(tvec)

fxerrmat = zeros(length(tvec), nr)
fyerrmat = zeros(length(tvec), nr)

for i in 1:nr
    fxerrmat[:,i] .= calc_normalized_error(aerostates.Fx[:,i][idxs], fxmat[:,i])
    fyerrmat[:,i] .= calc_normalized_error(aerostates.Fy[:,i][idxs], -fymat[:,i])
end

fxerravg = reshape(mean(abs.(fxerrmat[testidxs, :]), dims=1), :)
fyerravg = reshape(mean(abs.(fyerrmat[testidxs, :]), dims=1), :)

ferrplt = plot(xaxis="Radius (m)", yaxis="Time-Averaged Loading Error (%)", legend=:best)
plot!(ferrplt, rvec, fxerravg, lab=L"$F_x$", lw=2, seriescolor=1)
plot!(ferrplt, rvec, fyerravg, lab=L"$F_y$", lw=2, seriescolor=2)
display(ferrplt)




##### Plot 5 seconds of tip deflection error
t55idx = findfirst(x -> x==50.0, tvec[testidxs])
t60idx = findfirst(x -> x==60.0, tvec[testidxs])
tdx5 = tdxerr_norm[t55idx:t60idx]
tdy5 = tdyerr_norm[t55idx:t60idx]
tdz5 = tdzerr_norm[t55idx:t60idx]

t5 = tvec[testidxs][t55idx:t60idx]

td5plt = plot(xaxis="Time (s)", yaxis="Tip Deflection Error (%)")
plot!(td5plt, t5, tdx5, lab=L"$\delta x$", lw=2)
plot!(td5plt, t5, tdy5, lab=L"$\delta y$", lw=2)
plot!(td5plt, t5, tdz5, lab=L"$\delta z$", lw=2)
display(td5plt) #Doesn't profer any insight on the frequency behavior. 


errplt = plot(ferrplt, td5plt, layout=(1, 2), size=(1000, 400), leftmargin=5mm, bottommargin=5mm)
display(errplt)


nodeidx = 12

if nodeidx<10
    ofidx = "00$nodeidx"
elseif nodeidx<100
    ofidx = "0$nodeidx"
else
    ofidx = "$nodeidx"
end

alphaplt = plot(xaxis="Time (s)", yaxis="Angle of Attack (deg)", title="Node $nodeidx", legend=:best)
plot!(tvec, aerostates.alpha[:,nodeidx].*(180/pi), lab="R")
plot!(tvec, outs["AB1N"*ofidx*"Alpha"], lab="OF")
# display(alphaplt)


Vxind = outs["AB1N"*ofidx*"Vindx"]
Vyind = outs["AB1N"*ofidx*"Vindy"]
Vx = outs["AB1N"*ofidx*"Vx"]
Vy = outs["AB1N"*ofidx*"Vy"]
Wof = @. sqrt((Vx+Vxind)^2 + (Vy + Vyind)^2)

Uplt = plot(xaxis = "Time (s)", yaxis="Inflow Velocity (m/s)", title="Node $nodeidx", legend=:best)
plot!(tvec, Wof, lab="OF", linestyle=:dash, seriescolor=1)
plot!(tvec, aerostates.W[:,nodeidx], lab="R", seriescolor=2, alpha=0.7)
# display(Uplt)

phiplt = plot(xaxis="Time (s)", yaxis="Inflow Angle (deg)", title="Node $nodeidx", leg=:topleft)
plot!(tvec, aerostates.phi[:,nodeidx].*(180/pi), lab="R")
plot!(tvec, outs["AB1N"*ofidx*"Phi"], lab="OF")
# display(phiplt)


twist = @. (aerostates.phi[:,nodeidx] + -aerostates.alpha[:,nodeidx])*180/pi
twistof = @. outs["AB1N"*ofidx*"Phi"] + -outs["AB1N"*ofidx*"Alpha"]

twistplt = plot(xaxis="Time (s)", yaxis="Twist Angle (deg)", title="Node $nodeidx")
plot!(tvec, twist, lab="R")
plot!(tvec, outs["AB1N"*ofidx*"Theta"], lab="OF")
plot!(tvec, twistof, lab="OF", linestyle=:dash)
# display(twistplt)



aziplt = plot(tvec, azimuth, lab="AeroDyn", xaxis="Time (s)", yaxis="Azimuthal angle (deg)", leg=:top)
plot!(tvec, aerostates.azimuth.*(180/pi), lab="WATT")
plot!(tvec, outs["Azimuth"], lab="ElastoDyn", linestyle=:dash)
# display(aziplt)

Mx_of = zeros(nt)
Mx_r = zeros(nt)

for i = 1:nt
    Mx_of[i] = of.root_bending_moment(rvec, fxmat[i,:])
    Mx_r[i] = of.root_bending_moment(rvec, aerostates.Fx[idxs[i],:])
end





Mxerr = calc_error(Mx_r, Mx_of)[testidxs]
Mxerr_rel = calc_relative_error(Mx_r, Mx_of)[testidxs]
Mxerr_norm = calc_normalized_error(Mx_r, Mx_of)[testidxs]

println("Root Bending Moment Error Statistics (%): ")
@show maximum(abs.(Mxerr_rel)), mean(abs.(Mxerr_rel)), std(abs.(Mxerr_rel))
println("Root Bending Moment Absolute Error Statistics (N·m): ")
@show maximum(abs.(Mxerr)), mean(abs.(Mxerr)), std(abs.(Mxerr))
println("Root Bending Moment Normalized Error Statistics: ")
@show maximum(abs.(Mxerr_norm)), mean(abs.(Mxerr_norm)), std(abs.(Mxerr_norm))


Merrplt = plot(xaxis="Time (s)", yaxis="Root Bending Moment Error (%)", leg=false)
plot!(tvec[testidxs], Mxerr)
# display(Merrplt)









tiploads = plot(xaxis="Time (s)", yaxis="Tip Load (N)", legend=(0.2, -0.2), legend_columns=2) #
plot!(tiploads, tvec[idxs], fxmat[idxs,end], lab=L"$F_x$ - OpenFAST", seriescolor=1)
plot!(tiploads, tvec, aerostates.Fx[:,end], lab=L"$F_x$ - WATT.jl", seriescolor=4)
plot!(tiploads, tvec[idxs], -fymat[idxs,end], lab=L"$F_y$ - OpenFAST", seriescolor=2)
plot!(tiploads, tvec, aerostates.Fy[:,end], lab=L"$F_y$ - WATT.jl", seriescolor=5)
display(tiploads)


tipdefs = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)", legend=(0.1, -0.2), legend_columns=2) #
plot!(tipdefs, tvec, tipdx, lab=L"$\delta x$ - OpenFAST", seriescolor=1)
plot!(tipdefs, tvec, -tipdef_z, lab=L"$\delta x$ - WATT.jl", seriescolor=4)
plot!(tipdefs, tvec, tipdy, lab=L"$\delta y$ - OpenFAST", seriescolor=2)
plot!(tipdefs, tvec, tipdef_y, lab=L"$\delta y$ - WATT.jl", seriescolor=5)
plot!(tipdefs, tvec, tipdz, lab=L"$\delta z$ - OpenFAST", seriescolor=3)
plot!(tipdefs, tvec, tipdef_x, lab=L"$\delta z$ - WATT.jl", seriescolor=6)
display(tipdefs)


Mplt = plot(xaxis="Time (s)", leg=:best)
ylabel!(L"\parbox{15em}{Flapwise Root\\Bending Moment (N·m)}")
plot!(tvec, Mx_of, lab="OpenFAST", seriescolor=1)
plot!(tvec, Mx_r, lab="WATT.jl", seriescolor=4)
display(Mplt)






nothing