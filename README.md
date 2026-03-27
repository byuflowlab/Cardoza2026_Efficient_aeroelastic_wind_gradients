# Cardoza2025_Efficient_aeroelastic_wind_gradients
Scripts to reproduce results in Cardoza 2025 gradients for aeroelastic optimization paper. 

In order to run any of the scripts in this repository you must first install `DynamicStallModels`, `GXBeam`, `CCBlade`, `WATT`, and other packages. `DynamicStallModels` and `WATT` are not yet registered packages, so you may need to load them via url. The other packages you should be able to [load by activating the local environment](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project-1) (the manifest file is provided) and instantiating you should be able to get the other packages at the exact version we used them. 

```julia
add https://github.com/byuflowlab/DynamicStallModels.jl.git
add https://github.com/byuflowlab/WATT.jl.git
```

The scripts provided here are:
- verify.jl - A comparison of OpenFAST and WATT for a turbulent inflow case
- optimize.jl - An optimization of a simplified version of the NREL 5MW using sparse parallelized forward mode derivatives.
- Ndesignvars.jl - A comparison of the different forward differentiation techniques while varying the number of dense columns (chord and twist). 
- Nwind.jl - A comparison of the different forward differentiation techniques while varying the number of sparse columns (number of wind speeds). 
- *_reverse.jl - timing the reverse mode approach for Ndesignvars and Nwind. 
- random_start.jl - Run an optimization from a randomly perturbed design from the initial (designed to be run in a slurm job array). 
- parse_random_start_optimizations.jl - From a set of files labeled "perturbation_x.log" (The outputs of a random start optimization piped to a log file) analyze converged and failed designs. 

In order to run the all of the scripts, a module is provided in this repository that contains functions for differentiation and objective/constraint function analysis. If instantiating the environment above doesn't make the code accessible, simply `add ./UnsteadyOpt` in the Julia package manager when the local environment is active. Note that the optimization scripts require a license for SNOPT, but using SNOW you can use IPOPT (but you won't get the same results). 