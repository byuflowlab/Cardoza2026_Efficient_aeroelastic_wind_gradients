#=
Plot the unsteady deflections from the optimization code.

2/26/26 Adam Cardoza
=#

using Plots, DelimitedFiles


tvec = collect(0:0.05:100.0) #time vector, seconds (for fatigue analysis)

u = readdlm("deflections.csv", Float64) #Read the deflections from the CSV file
u = reshape(u, length(tvec), :, 3) #Reshape to (ntime, nelements, 3)


plt = plot(tvec, u[:, 36, 1], label="u_x", xlabel="Time (s)", ylabel="Deflection (m)", title="Unsteady Deflection of Element 36")
plot!(tvec, u[:, 36, 2], label="u_y")
plot!(tvec, u[:, 36, 3], label="u_z")
display(plt)

#This looks right 2/27/26

nothing
