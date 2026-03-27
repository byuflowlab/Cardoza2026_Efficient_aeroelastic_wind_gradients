# Cardoza2026_Efficient_aeroelastic_wind_gradients
Scripts to reproduce results in Cardoza 2026 paper ["Efficient derivative computation for unsteady fatigue-constrained nonlinear aero-structural wind turbine blade optimization"](https://doi.org/10.5194/wes-2026-10). 

In order to run any of the scripts in this repository you must first install Julia 1.11.x and the required packages. The `Manifest.toml` pins all dependencies to the exact versions used in the paper. To reproduce the environment:

```julia
# In the Julia REPL from the repo root:
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will download all packages at the exact versions we used, including the unregistered packages (`DynamicStallModels`, `GXBeam`, `WATT`, `GXBeamCS`, `OpenFASTTools`) which are pinned by git commit hash in the manifest.

**Note on Snopt:** The optimization scripts (`optimize.jl`, `random_start.jl`) use Snopt, which requires a commercial license for the underlying solver. If you do not have a Snopt license, you can use IPOPT via SNOW instead (results will differ). The `Snopt.jl` package itself will still be installed by `Pkg.instantiate()` — only running it requires a license.

The scripts provided here are:
- verify.jl - A comparison of OpenFAST and WATT for a turbulent inflow case
- optimize.jl - An optimization of a simplified version of the NREL 5MW using sparse parallelized forward mode derivatives.
- Ndesignvars.jl - A comparison of the different forward differentiation techniques while varying the number of dense columns (chord and twist). 
- Nwind.jl - A comparison of the different forward differentiation techniques while varying the number of sparse columns (number of wind speeds). 
- *_reverse.jl - timing the reverse mode approach for Ndesignvars and Nwind. 
- random_start.jl - Run an optimization from a randomly perturbed design from the initial (designed to be run in a slurm job array). 
- parse_random_start_optimizations.jl - From a set of files labeled "perturbation_x.log" (The outputs of a random start optimization piped to a log file) analyze converged and failed designs. 

In order to run the all of the scripts, a module is provided in this repository that contains functions for differentiation and objective/constraint function analysis. If instantiating the environment above doesn't make the code accessible, simply `add ./UnsteadyOpt` in the Julia package manager when the local environment is active.