module UnsteadyOpt

    using YAML, DelimitedFiles
    using GXBeam, GXBeamCS, CCBlade, OpenFASTTools
    import WATT
    using Base.Threads, StaticArrays, LinearAlgebra, FLOWMath, Statistics
    using RecipesBase

    using DiffResults, PolyesterForwardDiff, ForwardDiff, ReverseDiff
    using SparseArrays, SparseDiffTools

    import SpecialFunctions: gamma
    import GXBeam: nondual_value

    of = OpenFASTTools


    include("derivatives.jl") #How to compute derivatives for objectives and constraints
    include("utils.jl")
    include("aerodynamics.jl")
    include("composites.jl")
    include("structures.jl") 


    ### WISDEM cost model
    include("cost_table.jl")
    include("PPI.jl")
    include("costmodel.jl")

    include("parse.jl") #Parse SNOPT files for random start optimizations. 

end # module UnsteadyOpt
