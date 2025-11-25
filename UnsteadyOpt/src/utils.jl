
function designvar_name(idx, nx_cp, nx_cp_twist, nf_cp, nx_cp_w1, nx_cp_w2, nwind, nx)
    if idx == 1
        return "Movable first control point (r1)", 1
    elseif 2 <= idx <= 1+nx_cp
        return "Chord control point", idx-1
    elseif 1+nx_cp < idx <= 1+nx_cp+nx_cp_twist
        return "Twist control point", idx-(1+nx_cp)
    elseif 1+nx_cp+nx_cp_twist < idx <= 1+nx_cp+nx_cp_twist+5*nf_cp
        return "Segment scaling factor (f_segs)", idx-(1+nx_cp+nx_cp_twist)
    elseif 1+nx_cp+nx_cp_twist+5*nf_cp < idx <= 1+nx_cp+nx_cp_twist+5*nf_cp+nx_cp_w1
        return "Web 1 scaling factor", idx-(1+nx_cp+nx_cp_twist+5*nf_cp)
    elseif 1+nx_cp+nx_cp_twist+5*nf_cp+nx_cp_w1 < idx <= 1+nx_cp+nx_cp_twist+5*nf_cp+nx_cp_w1+nx_cp_w2
        return "Web 2 scaling factor", idx-(1+nx_cp+nx_cp_twist+5*nf_cp+nx_cp_w1)
    elseif 1+nx_cp+nx_cp_twist+5*nf_cp+nx_cp_w1+nx_cp_w2 < idx <= 1+nx_cp+nx_cp_twist+5*nf_cp+nx_cp_w1+nx_cp_w2+nwind
        return "Pitch", idx-(1+nx_cp+nx_cp_twist+5*nf_cp+nx_cp_w1+nx_cp_w2)
    elseif idx == nx
        return "TSR", 1
    else
        return "Unknown", idx
    end
end

"""
    get_constraint_name(idx, constraint_names, constraint_ends)

Return the constraint name corresponding to the given index `idx`.
"""
function get_constraint_name(idx, constraint_names, constraint_ends)
    for (i, end_idx) in enumerate(constraint_ends)
        start_idx = i == 1 ? 1 : constraint_ends[i-1] + 1
        if idx >= start_idx && idx <= end_idx
            return constraint_names[i]
        end
    end
    return "Unknown"
end


function get_near_violations(g, lg, ug, near_tol, constraint_names, constraint_ends)

    near_violations = Dict()
    for (i, cname) in enumerate(constraint_names)
        start_idx = i == 1 ? 1 : constraint_ends[i-1] + 1
        end_idx = constraint_ends[i]
        vals = g[start_idx:end_idx]
        lvals = lg[start_idx:end_idx]
        uvals = ug[start_idx:end_idx]
        near_idxs = findall(j -> (vals[j] < lvals[j] + near_tol) ||
                                (uvals[j] - near_tol < vals[j]), 1:length(vals))
        near_violations[cname] = near_idxs .+ (start_idx - 1)
    end
    return near_violations
end

# Given a constraint index within the buckling range, determine which element it occurs within
function buckling_constraint_to_element(idx, num_elements)
    # num_elements: vector of number of constraints per element (length = number of elements)
    # idx: global constraint index (1-based, within buckling range)
    cum_elements = cumsum(num_elements)
    elem = searchsortedfirst(cum_elements, idx)
    return elem
end

# Given a constraint index within the fatigue range, determine which cross-section it occurs within
function fatigue_constraint_to_cross_section(idx, fat_idxs, num_elements, Ncon)

    shifted_idx = idx - Ncon # Adjust index to be within fatigue constraints
    
    fat_num_elements = num_elements[fat_idxs]
    cum_fat = cumsum(fat_num_elements)
    elem_idx = searchsortedfirst(cum_fat, shifted_idx)
    cs_number = fat_idxs[elem_idx]
    element_number = shifted_idx - (elem_idx == 1 ? 0 : cum_fat[elem_idx-1])
    
    return cs_number, element_number
end



function get_dynamic_tip_deflection_constraints(g, N_bulk_full, ntwist, ntimecon)
    # Returns the dynamic tip deflection constraint values from g
    start_idx = N_bulk_full + 2 + ntwist + 1
    end_idx = start_idx + ntimecon - 1
    return g[start_idx:end_idx]
end

# # Check that buck_idxs is a list of continuous indices
function check_continuous_indices(buck_idxs)
    is_continuous = true
    for idx in 2:length(buck_idxs)
        idx_range = buck_idxs[idx]
        
        prev_range = buck_idxs[idx-1]
        if first(idx_range) != last(prev_range) + 1
            @warn("$(last(prev_range)) -> $(first(idx_range))")
            is_continuous = false
        end
    end
    return is_continuous
end

function print_constraint_stats(name, start_idx, end_idx, g)
    vals = g[start_idx:end_idx]
    minval = minimum(vals)
    minmag = vals[argmin(abs.(vals))]
    maxval = maximum(vals)
    meanval = mean(vals)
    medianval = median(vals)
    println("  $name:")
    println("    min:    $minval")
    println("    minmag: $minmag")
    println("    max:    $maxval")
    println("    mean:   $meanval")
    println("    median: $medianval")
end



#### Optimization functions
function get_list_property(list)
    return [list[i][j] for i in 1:length(list) for j in 1:length(list[i])]
end

function create_list_property(p)
    return [SMatrix{6,6}(p[i:i+35]) for i in 1:36:length(p)]
end

function check_constraints(g, lg, ug) #Todo: Test me (and below)
    gcheck = zeros(length(g))
    for i = eachindex(g)
        if g[i] < lg[i]
            gcheck[i] = g[i] - lg[i] #This should return a negative number if 
        elseif g[i] > ug[i]
            gcheck[i] = g[i] - ug[i]
        end
    end
    return gcheck
end



#### Utility functions

function write_dual(f, dual)
    s = "$(ForwardDiff.value.(dual)) \n"
    s = replace(s, "[" => "")
    s = replace(s, "]" => "")
    write(f, s)
end

function sweep_design_var(fun, x, idx, p, n, ng;atol=1e-6, verbose::Bool=false)

    if isapprox(x[idx], 0.0; atol)
        xu = atol
        xl = -atol
    else
        xu = x[idx]*(1+p)
        xl = x[idx]*(1-p)
    end

    xrange = range(xl, stop=xu, length=n)
    fvals = zeros(n)
    gvals = zeros(n, ng)

    # Threads.@threads for i in eachindex(xrange)
    for i in eachindex(xrange)
        verbose ? println("Evaluating design variable at index $i") : nothing
        xnew = deepcopy(x)
        xnew[idx] = xrange[i]
        gi = view(gvals, i, :)
        fvals[i] = fun(gi, xnew)
    end

    return fvals, gvals, xrange
end

function sweep_design_var(obj, con, x, idx, p, n, ng;atol=1e-6, verbose::Bool=false)

    if isapprox(x[idx], 0.0; atol)
        xu = atol
        xl = -atol
    else
        xu = x[idx]*(1+p)
        xl = x[idx]*(1-p)
    end

    xrange = range(xl, stop=xu, length=n)
    gvals = zeros(n, ng)
    fvals = zeros(n)
    Threads.@threads for i in eachindex(xrange)
    # for i in eachindex(xrange)
        verbose ? println("Evaluating design variable at index $i") : nothing
        xnew = deepcopy(x)
        xnew[idx] = xrange[i]
        gi = view(gvals, i, :)
        fvals[i] = obj(xnew)
        con(gi, xnew)
    end

    return fvals, gvals, xrange
end

function check_derivative(fvec, xrange, tol)
    dx = diff(xrange)
    dfdx = diff(fvec) ./ dx

    n = length(dfdx) - 1

    idxs = Int[]
    for i in 2:n
        abs(dfdx[i-1] - dfdx[i]) > tol * abs(dfdx[i]) ? push!(idxs, i) : nothing
        abs(dfdx[i+1] - dfdx[i]) > tol * abs(dfdx[i]) ? push!(idxs, i+1) : nothing
    end

    return unique(idxs)
end




function nearest_to(x, y)
    return argmin(abs.(x .- y))
end

function topN(A, N)
    inds = partialsortperm(vec(A), 1:N; rev=true)
    vals = vec(A)[inds]
    subs = CartesianIndices(A)[inds]
    return vals, subs
end

function continuous_index(vector)
    for i = 1:length(vector)-1
        if vector[i]+1 != vector[i+1]
            return false
        end
    end
    return true
end

function create_design_vector(x0, c, i)
    x = deepcopy(x0)
    x[i] = c
    return x
end



function get_seg_thickness(obj::ObjectiveFunction, chords::Vector{Float64}, twists::Vector{Float64}, fvec::Vector{Float64}, fw1::Vector{Float64}, fw2::Vector{Float64})
    ## Scale the layer thickness of the different sectors
    materials_i = GXBeamCS.Material{Float64}.(obj.composites.materials) # Make appropriate materials vector. 
    segs_webs = scale_segments(obj.composites.stations, fvec, fw1, fw2, materials_i, obj.idxs.num_segments)

    
    ### Get the cross sectional properties. 
    clt_list = [get_clt(
        obj.composites.stations[i].xaf, 
        obj.composites.stations[i].yaf, 
        chords[i], 
        0.0, 
        obj.distro.le_loc[i], 
        obj.composites.stations[i].xbreak, 
        obj.composites.stations[i].webloc, 
        segs_webs[i].segments, 
        segs_webs[i].webs
    ) for i in 1:obj.aero.nr]

    seg_thickness = zeros(Float64, obj.aero.nr, 7)

    for i in 1:5
        seg_thickness[:, i] = [GXBeamCS.get_section_thickness(clt_list[j], i) for j in 1:obj.aero.nr]
    end

    seg_thickness[7:end, 6] = [GXBeamCS.get_section_thickness(clt_list[j], 11) for j in 7:obj.aero.nr] #Todo: Not all of the cross sections have a 11th section, so this might not work for all cases.
    seg_thickness[7:end, 7] = [GXBeamCS.get_section_thickness(clt_list[j], 12) for j in 7:obj.aero.nr]

    return seg_thickness
end

function get_seg_thickness(obj::ObjectiveFunction, chords::Vector{Float64}, twists::Vector{Float64}, fvec::Vector{Float64})
    ## Scale the layer thickness of the different sectors
    materials_i = GXBeamCS.Material{Float64}.(obj.composites.materials) # Make appropriate materials vector. 
    segs_webs = scale_segments(obj.composites.stations, fvec, materials_i, obj.idxs.num_segments)

    
    ### Get the cross sectional properties. 
    clt_list = [get_clt(
        obj.composites.stations[i].xaf, 
        obj.composites.stations[i].yaf, 
        chords[i], 
        0.0, 
        obj.distro.le_loc[i], 
        obj.composites.stations[i].xbreak, 
        obj.composites.stations[i].webloc, 
        segs_webs[i].segments, 
        segs_webs[i].webs
    ) for i in 1:obj.aero.nr]

    seg_thickness = zeros(Float64, obj.aero.nr, 7)

    for i in 1:5
        seg_thickness[:, i] = [GXBeamCS.get_section_thickness(clt_list[j], i) for j in 1:obj.aero.nr]
    end

    seg_thickness[7:end, 6] = [GXBeamCS.get_section_thickness(clt_list[j], 11) for j in 7:obj.aero.nr] #Todo: Not all of the cross sections have a 11th section, so this might not work for all cases.
    seg_thickness[7:end, 7] = [GXBeamCS.get_section_thickness(clt_list[j], 12) for j in 7:obj.aero.nr]

    return seg_thickness
end

"""
    average_every_n(x, n)

Returns a vector containing the average of every `n` consecutive elements in `x`.
If the length of `x` is not a multiple of `n`, the last group will contain fewer elements.

**Arguments**
- x: Vector of numbers.
- n: Integer, number of points to average.

**Returns**
- Vector of averaged values.
"""
function average_every_n(x::AbstractVector, n::Int)
    result = similar(x, ceil(Int, length(x)/n))
    for i in 1:n:length(x)
        last_idx = min(i + n - 1, length(x))
        result[ceil(Int, i/n)] = mean(x[i:last_idx])
    end
    return result
end

export plot_blade_planform

@recipe function plot_blade_planform(objfun::ObjectiveFunction, chords)
    
    distro = objfun.distro
    rvec = distro.rvec
    leloc = distro.le_loc
    te_loc = @. 1 - leloc

    topy = chords.*leloc
    bottomy = chords.*te_loc

    # @show bottomy
    
    mainline = zeros(length(rvec))

    ribbon --> (bottomy, topy)
    linealpha --> 0.0

    rvec, mainline
end

export plot_cross_segs

@recipe function plot_cross_segs(obj::ObjectiveFunction, chords::Vector{Float64}, twists::Vector{Float64}, fvec::Vector{Float64}, fw1::Vector{Float64}, fw2::Vector{Float64}, idx::Int)

    ## Scale the layer thickness of the different sectors
    materials_i = GXBeamCS.Material{Float64}.(obj.composites.materials) # Make appropriate materials vector. 
    segs_webs = scale_segments(obj.composites.stations, fvec, fw1, fw2, materials_i, obj.idxs.num_segments)

    
    ### Get the cross sectional properties. 
    clt_list = [get_clt(
        obj.composites.stations[i].xaf, 
        obj.composites.stations[i].yaf, 
        chords[i], 
        0.0, 
        obj.distro.le_loc[i], 
        obj.composites.stations[i].xbreak, 
        obj.composites.stations[i].webloc, 
        segs_webs[i].segments, 
        segs_webs[i].webs
    ) for i in 1:obj.aero.nr]

    clt = clt_list[idx]

    # Find the max and min x and y coordinates for the cross section
    # xmin = 0.0
    # xmax = 0.0
    ymin = 0.0
    ymax = 0.0
    for sec in clt.sections
        # xmin = min(xmin, minimum(sec.y))
        # xmax = max(xmax, maximum(sec.y))
        ymin = min(ymin, minimum(sec.z))
        ymax = max(ymax, maximum(sec.z))
    end

    # @show length(clt.sections)

    # for i in eachindex(clt.sections)
    for i in 1:5
        @series begin

            xmin = minimum(clt.sections[i].y)
            xmax = maximum(clt.sections[i].y)

            x1 = [xmin, xmax]
            yl = -ones(2).*ymin
            yu = ones(2).*ymax

            # if i < 4
            #     seriescolor --> i
            # else
                seriescolor --> i+2
            # end

            ribbon --> (yl, yu)
            linealpha --> 0.0

            x1, zero(x1)
        end
    end

    # for i in 1:5
    # end

    # @show xmin, xmax, ymin, ymax

    counter = 1 #Element counter
    for i in eachindex(clt.sections)
        sec = clt.sections[i]
        ns = length(sec.y)

        for j in 1:ns-1 #Iterate over the number of elements in the section
            x = sec.y[j+1] - sec.y[j] #The x distance between the element end points
            y = sec.z[j+1] - sec.z[j] #The y distance between the element end points
            L = sqrt(x^2 + y^2)


            #Normal vector to the element
            nx = -y/L
            ny = x/L

            # @show nx, ny

            T = 0.0 #The total thickness travelled so far. 

            for k in eachindex(sec.laminate) #Iterate over the laminates in the section
                t = sec.laminate[k].t #The thickness of the laminate
                
                xp = [sec.y[j] + nx*T, sec.y[j+1] + nx*T, sec.y[j+1] + nx*(T+t), sec.y[j] + nx*(T+t), sec.y[j] + nx*T]
                yp = [sec.z[j] + ny*T, sec.z[j+1] + ny*T, sec.z[j+1] + ny*(T+t), sec.z[j] + ny*(T+t), sec.z[j] + ny*T]


                @series begin
                    label --> false
                    # if in(counter, highlight_counters)
                    #     seriescolor --> highlight
                    #     markershape --> :circle
            
                    # else
                    #     seriescolor --> :black
                    # end
                    linecolor --> :black

                    xp, yp
                end
                
                counter += 1 #Increment the element counter
                T += t #Increment the thickness
            end #End looping over laminates
        end #end looping over segment elements
    end #End looping over sections
end


@recipe function plot_cross_segs2(obj::ObjectiveFunction, chords::Vector{Float64}, twists::Vector{Float64}, fvec::Vector{Float64}, idx::Int)

    ## Scale the layer thickness of the different sectors
    materials_i = GXBeamCS.Material{Float64}.(obj.composites.materials) # Make appropriate materials vector. 
    segs_webs = scale_segments(obj.composites.stations, fvec, materials_i, obj.idxs.num_segments)

    
    ### Get the cross sectional properties. 
    clt_list = [get_clt(
        obj.composites.stations[i].xaf, 
        obj.composites.stations[i].yaf, 
        chords[i], 
        0.0, 
        obj.distro.le_loc[i], 
        obj.composites.stations[i].xbreak, 
        obj.composites.stations[i].webloc, 
        segs_webs[i].segments, 
        segs_webs[i].webs
    ) for i in 1:obj.aero.nr]

    clt = clt_list[idx]

    # Find the max and min x and y coordinates for the cross section
    # xmin = 0.0
    # xmax = 0.0
    ymin = 0.0
    ymax = 0.0
    for sec in clt.sections
        # xmin = min(xmin, minimum(sec.y))
        # xmax = max(xmax, maximum(sec.y))
        ymin = min(ymin, minimum(sec.z))
        ymax = max(ymax, maximum(sec.z))
    end

    # @show length(clt.sections)

    # for i in eachindex(clt.sections)
    for i in 1:5
        @series begin

            xmin = minimum(clt.sections[i].y)
            xmax = maximum(clt.sections[i].y)

            x1 = [xmin, xmax]
            yl = -ones(2).*ymin
            yu = ones(2).*ymax

            # if i < 4
            #     seriescolor --> i
            # else
                seriescolor --> i+2
            # end

            ribbon --> (yl, yu)
            linealpha --> 0.0

            x1, zero(x1)
        end
    end

    # for i in 1:5
    # end

    # @show xmin, xmax, ymin, ymax

    counter = 1 #Element counter
    for i in eachindex(clt.sections)
        sec = clt.sections[i]
        ns = length(sec.y)

        for j in 1:ns-1 #Iterate over the number of elements in the section
            x = sec.y[j+1] - sec.y[j] #The x distance between the element end points
            y = sec.z[j+1] - sec.z[j] #The y distance between the element end points
            L = sqrt(x^2 + y^2)


            #Normal vector to the element
            nx = -y/L
            ny = x/L

            # @show nx, ny

            T = 0.0 #The total thickness travelled so far. 

            for k in eachindex(sec.laminate) #Iterate over the laminates in the section
                t = sec.laminate[k].t #The thickness of the laminate
                
                xp = [sec.y[j] + nx*T, sec.y[j+1] + nx*T, sec.y[j+1] + nx*(T+t), sec.y[j] + nx*(T+t), sec.y[j] + nx*T]
                yp = [sec.z[j] + ny*T, sec.z[j+1] + ny*T, sec.z[j+1] + ny*(T+t), sec.z[j] + ny*(T+t), sec.z[j] + ny*T]


                @series begin
                    label --> false
                    # if in(counter, highlight_counters)
                    #     seriescolor --> highlight
                    #     markershape --> :circle
            
                    # else
                    #     seriescolor --> :black
                    # end
                    linecolor --> :black

                    xp, yp
                end
                
                counter += 1 #Increment the element counter
                T += t #Increment the thickness
            end #End looping over laminates
        end #end looping over segment elements
    end #End looping over sections
end

export myticks

function myticks(n, ytup...)
    y = vcat(ytup...)
    # @show minimum(y), maximum(y), n

    ticks = round.(Int, collect(range(minimum(y), stop=maximum(y), length=n))) 
    return ticks
end

"""
    cm2figsize(width_cm, height_cm)

Convert width and height in centimeters to figure size tuple for Plots.jl.
"""
function cm2figsize(width_cm, height_cm)
    px_per_inch = 96  # Plots.jl default DPI
    cm_per_inch = 2.54
    width_px = (width_cm / cm_per_inch) * px_per_inch
    height_px = (height_cm / cm_per_inch) * px_per_inch
    return ((round(Int, width_px)), round(Int, height_px))
end

"""
    in2figsize(width_in, height_in)

Convert width and height in inches to figure size tuple for Plots.jl.
"""
function in2figsize(width_in, height_in)
    px_per_inch = 96  # Plots.jl default DPI
    width_px = width_in * px_per_inch
    height_px = height_in * px_per_inch
    return (round(Int, width_px), round(Int, height_px))
end