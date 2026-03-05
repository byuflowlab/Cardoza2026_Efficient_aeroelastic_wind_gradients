"""
    ObjectiveFunction(distro, assembly, aero, structural, composites, idxs, scaling)

A structure for variable capture for the objective function. 

*Arguments*
- `distro::NamedTuple`: The precomp distribution.
- `assembly::GXBeam.Assembly`: The assembly information.
- `aero::NamedTuple`: The aerodynamic properties.
- `structural::NamedTuple`: The structural properties.
- `composites::NamedTuple`: The composite properties.
- `idxs::NamedTuple`: The indices for the design variables, constraints, and more.
- `scaling::NamedTuple`: The scaling factors for the design variables and constraints.
"""
struct ObjectiveFunction{T1, T2, T3, T4, T5, T6, T7}
    distro::T1
    assembly::T2
    aero::T3
    structural::T4
    composites::T5
    idxs::T6
    scaling::T7
end

struct ConstraintFunction{T1, T2, T3, T4, T5, T6, T7}
    distro::T1
    assembly::T2
    aero::T3
    structural::T4
    composites::T5
    idxs::T6
    scaling::T7
end

struct Deriv{T1, T2, T3}
    obj::T1
    y::T2
    dx::T3
end

"""
    Deriv(obj::ObjectiveFunction, sizes::Tuple{Int, Int})

Create a derivative object for the given objective function and sizes.

*Arguments*
- `obj::ObjectiveFunction`: The objective function for which derivatives are computed.
- `sizes::Tuple{Int, Int}`: A tuple containing the number of design variables and the number of constraints.
"""
function Deriv(obj::ObjectiveFunction, sizes::Tuple{Int, Int})

    nx = sizes[1]
    ng = sizes[2] + 1

    y = zeros(ng)
    dx = zeros(ng, nx)

    return Deriv(obj, y, dx)
end

function (deriv::Deriv)(g, df, dg, x)
    nx = length(x)
    ng = length(g)

    function combine!(y, x)
        gg = view(y, 2:ng+1)
        y[1] = deriv.obj(gg, x)
    end

    PolyesterForwardDiff.threaded_jacobian!(combine!, deriv.y, deriv.dx, x, ForwardDiff.Chunk(8)) #Todo: Make the chunk size a parameter. 



    df .= view(deriv.dx, 1, :)
    dg .= view(deriv.dx, 2:ng+1, :)
    g  .= view(deriv.y, 2:ng+1)
    return deriv.y[1]
end




struct SparseDeriv{T1, T2, T3, T4, TT, TR, TS, TC}
    obj::T1 #Objective object (functor)
    con::T2 #Constraints object (functor)
    y::T3 #gradient
    dx::T4 #todo: I think this is unused. 
    tape::TT #Objective tape
    result::TR #objective gradient DiffResults cache #todo: convert to DifferentiationInterface? 
    sparseM::TS #gradient backend
    colorvec::TC #gradient cache
    chunksize::Int
end

"""
    SparseDeriv(obj::ObjectiveFunction, sizes::Tuple{Int, Int})

Create a derivative object for the given objective function and sizes.

*Arguments*
- `obj::ObjectiveFunction`: The objective function for which derivatives are computed.
- `sizes::Tuple{Int, Int}`: A tuple containing the number of design variables and the number of constraints.
"""
function SparseDeriv(obj::ObjectiveFunction, con::ConstraintFunction, nx, ng, sparsity_pattern, chunksize, x0)

    y = zeros(ng) #constraints
    dx = zeros(ng, nx) #constraint jacobian

    # x0 = zeros(nx)

    ### Objective cache
    tape = ReverseDiff.GradientTape(obj, x0) 
    compiledtape = ReverseDiff.compile(tape)
    @warn "Using ReverseDiff for objective gradient. Some regions of the design space break ReverseDiff. Consider using SparseForwardDeriv for better performance."

    result = DiffResults.GradientResult(x0)

    ### Constraints cache
    colorvec = matrix_colors(sparsity_pattern)

    return SparseDeriv(obj, con, y, dx, compiledtape, result, sparsity_pattern, colorvec, chunksize) 
end

function (deriv::SparseDeriv)(g, df, dg, x)

    ### Objective function
    # println("computing gradient...")
    ReverseDiff.gradient!(deriv.result, deriv.tape, x) 
    # println("Fetching gradient...")
    df .= DiffResults.gradient(deriv.result)

    ### Constraint functions
    # println("Computing jacobian")
    dg .= colored_jacobian_polyesterforwarddiff(deriv.con, g, x, deriv.sparseM; colorvec=deriv.colorvec, chunksize=deriv.chunksize)

    return DiffResults.value(deriv.result)
end



struct SparseForwardDeriv{T1, T2, T3, T4, TT, TR, TS, TC}
    obj::T1 #Objective object (functor)
    con::T2 #Constraints object (functor)
    y::T3 #gradient
    dx::T4 #todo: I think this is unused. 
    tape::TT #Objective tape
    result::TR #objective gradient DiffResults cache #todo: convert to DifferentiationInterface? 
    sparseM::TS #gradient backend
    colorvec::TC #gradient cache
    chunksize::Int
end

"""
    SparseDeriv(obj::ObjectiveFunction, sizes::Tuple{Int, Int})

Create a derivative object for the given objective function and sizes.

*Arguments*
- `obj::ObjectiveFunction`: The objective function for which derivatives are computed.
- `sizes::Tuple{Int, Int}`: A tuple containing the number of design variables and the number of constraints.
"""
function SparseForwardDeriv(obj::ObjectiveFunction, con::ConstraintFunction, nx, ng, sparsity_pattern, chunksize, x0)

    y = zeros(ng) #constraints
    dx = zeros(ng, nx) #constraint jacobian

    # x0 = zeros(nx)

    ### Objective cache
    tape = ReverseDiff.GradientTape(obj, x0) 
    compiledtape = ReverseDiff.compile(tape)

    result = DiffResults.GradientResult(x0)

    ### Constraints cache
    colorvec = matrix_colors(sparsity_pattern)

    return SparseForwardDeriv(obj, con, y, dx, compiledtape, result, sparsity_pattern, colorvec, chunksize) 
end

function (deriv::SparseForwardDeriv)(g, df, dg, x)

    ### Objective function
    # # println("computing gradient...")
    # ReverseDiff.gradient!(deriv.result, deriv.tape, x) 
    # # println("Fetching gradient...")
    # df .= DiffResults.gradient(deriv.result)

    ForwardDiff.gradient!(deriv.result, deriv.obj, x)
    df .= DiffResults.gradient(deriv.result)

    ### Constraint functions
    # println("Computing jacobian")
    dg .= colored_jacobian_polyesterforwarddiff(deriv.con, g, x, deriv.sparseM; colorvec=deriv.colorvec, chunksize=deriv.chunksize)

    return DiffResults.value(deriv.result)
end

struct ForwardDeriv{T1, T2, T3, TR}
    obj::T1 #Objective object (functor)
    con::T2 #Constraints object (functor)
    y::T3 #gradient
    result::TR #objective gradient DiffResults cache 
    chunksize::Int
end

"""
    ForwardDeriv(obj::ObjectiveFunction, sizes::Tuple{Int, Int})

Create a derivative object for the given objective function and sizes.

*Arguments*
- `obj::ObjectiveFunction`: The objective function for which derivatives are computed.
- `sizes::Tuple{Int, Int}`: A tuple containing the number of design variables and the number of constraints.
"""
function ForwardDeriv(obj::ObjectiveFunction, con::ConstraintFunction, nx, ng, chunksize, x0)
    y = zeros(ng) #constraints

    ### Objective cache
    result = DiffResults.GradientResult(x0)


    return ForwardDeriv(obj, con, y, result, chunksize) 
end

function (deriv::ForwardDeriv)(g, df, dg, x)

    ### Objective function
    ForwardDiff.gradient!(deriv.result, deriv.obj, x)
    df .= DiffResults.gradient(deriv.result)

    ### Constraint functions
    PolyesterForwardDiff.threaded_jacobian!(deriv.con, g, dg, x, ForwardDiff.Chunk(deriv.chunksize))

    return DiffResults.value(deriv.result)
end

struct FiniteDeriv{T1, T2, T3, T4, T5}
    obj::T1
    con::T2
    y::T3
    dx::T4
    sp::T5
    colors::Vector{Int}
end

"""
    Deriv(obj::ObjectiveFunction, sizes::Tuple{Int, Int})

Create a derivative object for the given objective function and sizes.

*Arguments*
- `obj::ObjectiveFunction`: The objective function for which derivatives are computed.
- `sizes::Tuple{Int, Int}`: A tuple containing the number of design variables and the number of constraints.
"""
function FiniteDeriv(obj::ObjectiveFunction, con::ConstraintFunction, nx, ng, sparsity_pattern)

    y = zeros(ng)
    dx = zeros(ng, nx)

    colors = matrix_colors(sparsity_pattern)

    return FiniteDeriv(obj, con, y, dx, sparsity_pattern, colors)
end

function (deriv::FiniteDeriv)(g, df, dg, x)
    ng = length(g)

    f = deriv.obj(x)
    df .= center_gradient_parallel(deriv.obj, x)

    dg .= colored_parallel_finite_jacobian!(deriv.con, deriv.y, x, deriv.sp; colorvec=deriv.colors)

    # df .= view(deriv.dx, 1, :)
    # dg .= view(deriv.dx, 2:ng+1, :)

    return f
end

function relative_step(x, eps_rel)
    return eps_rel.*(max.(abs.(x), 1.0))
end



function center_jacobian_parallel(fun, x; eps_rel=1e-8, eps=relative_step(x, eps_rel), verbose::Bool=false)
    ni = length(x)
    y0 = fun(deepcopy(x))
    no = length(y0)

    J = zeros(no, ni)
    
    Threads.@threads for i = 1:ni
        verbose ? println("Column $i") : nothing
        xi = deepcopy(x)
        xi[i] += eps[i]
        yp = fun(xi)
        xi[i] -= 2*eps[i]
        ym = fun(xi)
        J[:, i] = @. (yp - ym) / (2*eps[i])
    end

    return J
end

function center_gradient_parallel(fun, x; eps_rel=1e-8, eps=relative_step(x, eps_rel), verbose::Bool=false)
    ni = length(x)

    df = zeros(ni)
    
    Threads.@threads for i = 1:ni
        verbose ? println("Column $i") : nothing
        xi = deepcopy(x)
        xi[i] += eps[i]
        yp = fun(xi)
        xi[i] -= 2*eps[i]
        ym = fun(xi)
        df[i] = (yp - ym) / (2*eps[i])
    end

    return df
end

function center_jacobian_parallel!(J, fun, x; eps_rel=1e-8, eps=relative_step(x, eps_rel), verbose::Bool=false)
    ni = length(x)
    
    Threads.@threads for i = 1:ni
        verbose ? println("Column $i") : nothing
        xi = deepcopy(x)
        xi[i] += eps[i]
        yp = fun(xi) #We don't want to pre-allocate these because we want them to be different memory for each thread.
        xi[i] -= 2*eps[i]
        ym = fun(xi)
        J[:, i] = @. (yp - ym) / (2*eps[i])
    end

    return J
end





"""
    map_inputs(xc, x, colors)
Map the input duals `xc` to the full input `x` using the `colors` vector.

**Arguments:**
- `xc`: A vector of dual numbers valued at one in the compressed space.
- `x`: The full input vector.
- `colors`: A vector indicating the color of each element in the sparse matrix.
"""
function map_inputs(xc, x, colorvec)
    # Create a new vector of duals with proper seeding
    xm = Vector{eltype(xc)}(undef, length(x))
    
    for i in eachindex(x)
        color = colorvec[i]
        xm[i] = typeof(xc[color])(x[i], ForwardDiff.partials(xc[color]))
    end
    
    return xm
end


function colored_jacobian_forwarddiff(f!, y, x, sparsity_pattern; colorvec=matrix_colors(sparsity_pattern))

    function wrap!(yy, xc)
        xr = map_inputs(xc, x, colorvec)
        f!(yy, xr)
        nothing
    end

    colors = unique(colorvec)

    xc = ones(length(colors)) #Can cache

    Jc = ForwardDiff.jacobian(wrap!, y, xc)

    # return Jc
    rows, cols, vals = SparseArrays.findnz(sparsity_pattern)
    ones_sparse = sparse(rows, cols, ones(length(vals)), size(sparsity_pattern, 1), size(sparsity_pattern, 2)) #Can move outside

    Jout = zeros(length(y), length(x))

    for k in eachindex(colors)
        mask = colorvec .== colors[k]
        backidxs = findall(mask) # get the indices of the mask

        ji = view(Jc, :, k)

        #Map the reduced Jacobian back to the full matrix
        for j in backidxs
            for i in eachindex(ji) 
                Jout[i, j] += ji[i]*ones_sparse[i, j] #values appears to be a tuple.
            end
        end
    end

    return Jout
end


struct WrappedFunction{F}
    f::F
    x::Vector{Float64}
    colorvec::Vector{Int}
end

function (wrap::WrappedFunction)(yy, xc)
    xr = map_inputs(xc, wrap.x, wrap.colorvec)
    wrap.f(yy, xr)
    nothing
end

function update_wrap!(wrap::WrappedFunction, new_x::Vector{Float64})
    wrap.x .= new_x
end

struct ColorCache{TS, TF, TJ, TC}
    sp::TS # Sparse matrix
    colorvec::Vector{Int} # Vector of colors for each column
    xc::Vector{Float64} # Cached values for the colors
    colors::Vector{Int} # Unique colors
    wrap::TF
    Jc::TJ
    config::TC
end

function prep_color_cache(sparsity_pattern, f!, y, x)
    colorvec = matrix_colors(sparsity_pattern)
    colors = unique(colorvec)
    xc = ones(length(colors)) 

    rows, cols, vals = SparseArrays.findnz(sparsity_pattern)
    ones_sparse = sparse(rows, cols, ones(length(vals)), size(sparsity_pattern, 1), size(sparsity_pattern, 2)) 

    wrap = WrappedFunction(f!, x, colorvec)

    cfg = ForwardDiff.JacobianConfig(wrap, y, xc)
    Jc = zeros(length(y), length(xc))

    return ColorCache(ones_sparse, colorvec, xc, colors, wrap, Jc, cfg)
end


function colored_jacobian_forwarddiff(cache::ColorCache, Jout, y, x)

    update_wrap!(cache.wrap, x) # Update the wrapped function with the new x values

    ForwardDiff.jacobian!(cache.Jc, cache.wrap, y, cache.xc, cache.config)

    for k in eachindex(cache.colors)
        mask = cache.colorvec .== cache.colors[k]
        backidxs = findall(mask) # get the indices of the mask

        ji = view(cache.Jc, :, k)

        #Map the reduced Jacobian back to the full matrix
        for j in backidxs
            for i in eachindex(ji) 
                Jout[i, j] += ji[i]*cache.sp[i, j] #values appears to be a tuple.
            end
        end
    end

    return Jout
end



function colored_jacobian_polyesterforwarddiff(f!, y, x, sparsity_pattern; colorvec=matrix_colors(sparsity_pattern), chunksize=8)

    function wrap!(yy, xc)
        xr = map_inputs(xc, x, colorvec)
        f!(yy, xr)
        nothing
    end

    colors = unique(colorvec)

    xc = ones(length(colors)) #Can cache
    Jc = zeros(length(y), length(colors))

    PolyesterForwardDiff.threaded_jacobian!(wrap!, y, Jc, xc, ForwardDiff.Chunk(chunksize))

    # return Jc
    rows, cols, vals = SparseArrays.findnz(sparsity_pattern)
    ones_sparse = sparse(rows, cols, ones(length(vals)), size(sparsity_pattern, 1), size(sparsity_pattern, 2)) #Can move outside

    Jout = zeros(length(y), length(x))

    for k in eachindex(colors) #Todo: This can put in a function and unified across methods. 
        mask = colorvec .== colors[k]
        backidxs = findall(mask) # get the indices of the mask

        ji = view(Jc, :, k)

        #Map the reduced Jacobian back to the full matrix
        for j in backidxs
            for i in eachindex(ji) 
                Jout[i, j] += ji[i]*ones_sparse[i, j] #values appears to be a tuple.
            end
        end
    end

    return Jout
end

function complex_gradient(f, x; h=1e-20)
    # Compute the complex step derivative
    g = zeros(length(x))
    xp = Complex.(x)

    for i in eachindex(x)
        xp[i] += im*h
        g[i] = imag(f(xp)) / h
        xp[i] = x[i]
    end
    return g
end

function complex_jacobian!(f, y, x; h=1e-20)
    # Compute the complex step derivative
    ni = length(x)
    no = length(y)

    J = zeros(no, ni)
    xp = Complex.(x)

    for i in eachindex(x)
        xp[i] += im*h
        f(y, xp)
        J[:, i] = @. imag(y) / h
        xp[i] = x[i]
    end

    return J
end




function colored_parallel_finite_jacobian!(f!, y, x, sparsity_pattern; colorvec=matrix_colors(sparsity_pattern), eps_rel=1e-8, eps=relative_step(x, eps_rel))

    ni = length(x)
    no = length(y)

    colors = unique(colorvec)

    rows, cols, vals = SparseArrays.findnz(sparsity_pattern)
    ones_sparse = sparse(rows, cols, ones(length(vals)), size(sparsity_pattern, 1), size(sparsity_pattern, 2)) #Can move outside

    Jout = zeros(no, ni)

    Threads.@threads for k in eachindex(colors)
        xi = deepcopy(x)
        # Efficiently find indices for the current color and perturb them
        mask = colorvec .== colors[k]
        backidxs = findall(mask) # get the indices of the mask
        @inbounds for i in backidxs
            xi[i] += eps[i]
        end
        yp = similar(y)
        f!(yp, xi) #We don't want to pre-allocate these because we want them to be different memory for each thread.

        @inbounds for i in backidxs
            xi[i] -= 2*eps[i]
        end
        ym = similar(y)
        f!(ym, xi)


        Ji = yp .- ym

        #Map the reduced Jacobian back to the full matrix
        for j in backidxs
            for i in eachindex(Ji) 
                Jout[i, j] += Ji[i]*ones_sparse[i, j]/(2*eps[j])
            end
        end
    end


    return Jout
end


"""
    compare_sparsity(A, B)
Compare the sparsity patterns of two matrices A and B.

**Arguments**
- `A`: The first sparse matrix to compare.
- `B`: The second sparse matrix to compare.

**Returns**
- `crA::Vector{Tuple{Int, Int}}`: The row and column indices of the non-zero elements in A, but not in B
- `crB::Vector{Tuple{Int, Int}}`: The row and column indices of the non-zero elements in B, but not in A.
"""
function compare_sparsity(A::SparseMatrixCSC, B::SparseMatrixCSC)

    rowA, colA, _ = SparseArrays.findnz(A)
    rowB, colB, _ = SparseArrays.findnz(B)

    rca = [(i, j) for (i, j) in zip(rowA, colA)]
    rcb = [(i, j) for (i, j) in zip(rowB, colB)]

    return setdiff(rca, rcb), setdiff(rcb, rca)
end



"""
    sample_space(lx, ux)

Generate a random vector `x` where each element is within its respective bounds.

*Arguments*
- `lx::AbstractVector`: Lower bounds for each variable.
- `ux::AbstractVector`: Upper bounds for each variable.

*Returns*
- `x::Vector{Float64}`: A vector with elements sampled uniformly within [lx[i], ux[i]].
"""
function sample_space(xold, p, lx::AbstractVector, ux::AbstractVector)
    @assert length(xold) == length(lx) == length(ux)

    xnew = similar(xold)
    for i in eachindex(lx)
        delta = (2*rand() - 1) * p * abs(xold[i]) # random value in [-p, p] * |xold[i]|
        # @show delta, xold[i], lx[i], ux[i]
        xnew[i] = xold[i] + delta

        if xnew[i] < lx[i]
            xnew[i] = lx[i]
        elseif xnew[i] > ux[i]
            xnew[i] = ux[i]
        end
    end

    return xnew
end



function scale_gradient(df; smax::Int=5, bottom::Float64=0.9, top::Float64=10.0)
    individual_scaling = ones(length(df))
    for i in eachindex(df)
        for j in 1:smax
            if abs(df[i]*individual_scaling[i])<bottom
                individual_scaling[i] *= 10
            elseif abs(df[i]*individual_scaling[i])>top
                individual_scaling[i] /= 10
            end
        end
    end
    return individual_scaling
end

function scale_jacobian(J; smax::Int=5, bottom::Float64=0.9, top::Float64=10.0)
    ng = size(J, 1)
    individual_scaling = ones(ng)
    for i in 1:ng
        for j in 1:smax
            if maximum(abs.(J[i, :]*individual_scaling[i]))<bottom
                individual_scaling[i] *= 10
            elseif maximum(abs.(J[i, :]*individual_scaling[i]))>top
                individual_scaling[i] /= 10
            end
        end
    end
    return individual_scaling
end


function get_vector_tracked_reals(xargs...)
    n::Int = 0
    for x in xargs
        n += length(x)
    end
    TF = promote_type(eltype.(xargs)...)
    xtracked = Vector{TF}(undef, n)
    i = 1
    for x in xargs
        for xi in x
            xtracked[i] = xi
            i += 1
        end
    end
    return xtracked
end

function get_tracked_array(example, xargs...)
    n::Int = 0
    for x in xargs
        n += length(x)
    end
    
    xtracked = similar(example, n)
    i = 1
    for x in xargs
        for xi in x
            xtracked[i] = xi
            i += 1
        end
    end
    return xtracked
end