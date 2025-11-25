#=
Parse Snopt files

=#



"""
"""
function clean_line(line, index)
    li =
        parse.(
            Float64,
            split(
                strip(replace(replace(line[index], "(" => " ", ")" => " "), "  " => " ")),
                " ",
            ),
        )

    return [li; line[(index + 1):(end - length(li) + 1)]]
end

"""
    parse_snopt_summary_file(filename)

Parse SNOPT summary file.

**Arguments:**
- `filename::String` : file to parse.

**Returns:**
- `snopt_summary::NamedTuple` : named tuple containing pertinent information.

**Contents of `snopt_summary`:**
- `:niterations::Int` : total number of iterations,
- `:nmajor::Int` : number of major iterations,
- `:function_calls::Int` : total number of function calls,
- `:major_iteration::Array{Int}` : major iteration number,
- `:minor_iterations::Array{Int}` : minor iterations in each major iteration,
- `:step_size::Array{Float}` : step size,
- `:nCons::Array{Int}` : nCons,
- `:feasibility::Array{Float}` : feasibility value for each iteration,
- `:optimality::Array{Float}` : optimality value for each iteration,
- `:merit_function::Array{Float}` : objective value for each iteration,
- `:nSs::Array{Int}` : nS value for each iteration (what is this?),
- `:exit::String` : exit message,
- `:info::String` : info message,
"""
function parse_snopt_summary_file(filename)

    #import summary file as set of lines
    # f = open(filename)
    lines = DelimitedFiles.readdlm(filename)
    # close(f)

    #initialize outputs
    minor_iterations = Int[]
    major_iteration_number = Int[]
    step_size = Float64[]
    nCons = Int[]
    feasibility = Float64[]
    optimality = Float64[]
    merit_function = Float64[]
    nSs = Int[]
    total_iterations = 0
    major_iterations = 0
    function_calls = 0
    exit_message = ""
    info_message = ""

    minor_flag = true
    #loop through lines
    for i in 2:length(lines[:, 1]) #ignore all the header stuff
        #check if it's one of the final lines that you need a number from
        if join(lines[i, 1:2], " ") == "SNOPTA EXIT"
            #save exit message
            exit_message = strip(join(lines[i, 3:end], " "))

        elseif join(lines[i, 1:2], " ") == "SNOPTA INFO"
            #save info message
            info_message = strip(join(lines[i, 3:end], " "))

        elseif join(lines[i, 1:3], " ") == "No. of iterations"
            #save total number of iterations
            total_iterations = lines[i, 4]

        elseif join(lines[i, 1:4], " ") == "No. of major iterations"
            #save number of major iterations
            major_iterations = lines[i, 5]

        elseif join(lines[i, 1:4], " ") == "User function calls (total)"
            #save number of function calls
            function_calls = lines[i, 5]

        elseif lines[i, 1] == "Major"
            minor_flag = false

        elseif lines[i, 1] == "Minor" || lines[i, 2] == 0
            minor_flag = true

        elseif typeof(lines[i, 1]) == Int
            if !minor_flag
                #modify vector if zeroth major iteration since step size is missing
                if lines[i, 1] == 0 || typeof(lines[i, 3]) != Float64
                    lines[i, 3:end] = [0.0; lines[i, 3:(end - 1)]]
                end

                #there's a string in here. need to clean out the parentheses.
                if typeof(lines[i, 5]) == SubString{String}
                    lines[i, 5:end] = clean_line(lines[i, :], 5)
                end
                if typeof(lines[i, 6]) == SubString{String}
                    lines[i, 6:end] = clean_line(lines[i, :], 6)
                end

                #don't save repeated stuff at the end
                if if (length(major_iteration_number) > 0)
                    (major_iteration_number[end] != lines[i, 1])
                else
                    true
                end

                    #push all the datalines
                    push!(major_iteration_number, lines[i, 1]) #major iteration number is first value

                    push!(minor_iterations, lines[i, 2]) #minor iterations is 2nd value

                    push!(step_size, lines[i, 3]) #step size is 3rd value

                    push!(nCons, lines[i, 4]) #nCons (what are these?) is 4th value

                    #check if feasibility is in parentheses (zero)
                    if typeof(lines[i, 5]) == Float64
                        push!(feasibility, lines[i, 5]) #feasibility is 5th value
                    else
                        push!(feasibility, 0.0)
                    end

                    push!(optimality, lines[i, 6]) #optimality is 6th value

                    push!(merit_function, lines[i, 7]) #objective value is 7th value

                    if typeof(lines[i, 8]) == SubString{String}
                        push!(nSs, 0)
                    elseif typeof(lines[i, 8]) == Float64
                        push!(nSs, 0)
                    else
                        push!(nSs, lines[i, 8]) #nS (what is this?) is 8th value
                    end
                end
            end
        end
    end

    snopt_summary = (
        niterations=total_iterations,
        nmajor=major_iterations,
        function_calls=function_calls,
        major_iteration=major_iteration_number,
        minor_iterations=minor_iterations,
        step_size=step_size,
        nCons=nCons,
        feasibility=feasibility,
        optimality=optimality,
        merit_function=merit_function,
        nSs=nSs,
        exit=exit_message,
        info=info_message,
    )

    return snopt_summary
end

function rmspaces(line)
    line = replace(line, "\t" => " ")
    # println("line: ", line)
    while occursin("  ", line)
        line = replace(line, "  "=>" ")
        # println("line: ", line)
    end

    return line
end

function parseentry(word; tryint=false)
    out = nothing
    if tryint
        out = tryparse(Int, word)
    end

    if isnothing(out)
        out = tryparse(Float64, word)
    end

    if isnothing(out)
        out = tryparse(Bool, word)
    end

    if isnothing(out)
        if contains(word, ",")
            try
                out = vec(readdlm(IOBuffer(word), ','))
            catch
                out = nothing
            end
        else
            out = nothing
        end
    end

    if isnothing(out)
        out = word
    end

    return out
end



function parseline(line_initial)
    line = removecomment(line_initial)
    # println("\"", line, "\"")
    entryend = findlast(" ", line)[1]-1 #Todo. This will cause problem if the requested line is a vector. -> I'll make a function that does this functionality. Inside that it'll see if this space occurs before or after a comma. -> Or just use findlast().... because I've used strip(), then the last space that occurs, should be just before the key. 
    entry = parseentry(line[1:entryend])
    key = nospaces(line[entryend+1:end])
    # println("\"", key, "\"")

    return key, entry
end

function parsefloat(line)
    if occursin("(", line)
        line = replace(line, "(" => "")
        line = replace(line, ")" => "")
        val = tryparse(Float64, line)
        return !isnothing(val) ? -val : 0
    else
        val = tryparse(Float64, line)
        return !isnothing(val) ? val : 0
    end
end


function parse_snopt_summary(filename)
    fi = open(filename, "r")
    lines = readlines(fi)
    close(fi)

    start_idx = 0
    nl = length(lines)
    for i in eachindex(lines)
        if occursin("Major Minor", lines[i])
            start_idx = i + 1
            break
        end
    end

    end_idx = 0
    for i = start_idx:nl
        if occursin("SNOPTA EXIT", lines[i])
            end_idx = i - 2
            while length(lines[end_idx]) < 10
                end_idx -= 1
            end
            break
        end
    end


    mylines = lines[start_idx:end_idx]

    
    idxs = Int[]
    i = 1
    # for i in eachindex(mylines)
    while i <= length(mylines)

        if mylines[i] == ""
            i += 1
        elseif occursin("Major Minor", mylines[i])
            i += 1
        elseif occursin("QP mult", mylines[i]) #Trying to get into feasible space
            skipsection = true
            while skipsection
                if occursin("Major Minor", mylines[i])
                    skipsection = false
                end
                i += 1
            end
        else
            @show mylines[i]
            push!(idxs, i)
            i += 1
        end

    end

    mylines = mylines[idxs]


    major = zeros(Int, length(mylines))
    minor = zeros(Int, length(mylines))
    step = zeros(length(mylines))
    ncons = zeros(Int, length(mylines))
    feasibility = zeros(length(mylines))
    optimality = zeros(length(mylines))
    merit = zeros(length(mylines))
    penalty = zeros(length(mylines))

    for i in eachindex(mylines)
        # @show i, mylines[i]
        line = mylines[i]
        # Major iter
        major[i] = parse(Int, line[1:6])

        # Minor iter
        minor[i] = parse(Int, line[8:14])

        # Step size
        step[i] = parsefloat(line[16:24])

        # nCons
        if isempty(strip(line[26:30]))
            ncons[i] = 0
        else
            ncons[i] = parse(Int, line[26:30])
        end
        # ncons[i] = parse(Int, line[26:30])

        # Feasibility
        feasibility[i] = parsefloat(line[31:39]) 

        # Optimality
        @show i, line[40:49]
        optimality[i] = parsefloat(line[40:49]) 

        # Merit function
        merit[i] = parsefloat(line[50:62]) 

        # Penalty
        if length(line) < 71
            penalty[i] = 0
            continue
        else
            penalty[i] = parsefloat(line[71:77]) 
        end
        
        
    end

    return (;major=major,
        minor=minor,
        step=step,
        ncons=ncons,
        feasibility=feasibility,
        optimality=optimality,
        merit=merit,
        penalty=penalty)
end


"""
Parse a perturbation log file to extract optimization results.

Returns a named tuple with:
- success: Bool indicating if optimization finished successfully
- x0: Initial design vector
- f0: Initial objective value  
- xopt: Optimal design vector
- fopt: Optimal objective value
- perturbation_num: Perturbation number from filename
- start_time: Optimization start timestamp
"""
function parse_perturbation_log(logfile_path::String)
    if !isfile(logfile_path)
        error("Log file not found: $logfile_path")
    end
    
    # Initialize return values
    success = false
    x0 = Float64[]
    f0 = NaN
    xopt = Float64[]
    fopt = NaN
    perturbation_num = 0
    start_time = ""
    
    # Extract perturbation number from filename
    filename = basename(logfile_path)
    if occursin(r"perturbation_(\d+)", filename)
        perturbation_num = parse(Int, match(r"perturbation_(\d+)", filename).captures[1])
    end
    
    lines = readlines(logfile_path)
    
    for (i, line) in enumerate(lines)
        # Extract start time and perturbation info
        if startswith(line, "running randomized_coe_opt.jl for perturbation")
            parts = split(line)
            if length(parts) >= 6
                start_time = parts[end]
            end
        end
        
        # Check for successful completion
        if occursin("Optimization complete.", line)
            success = true
        end

        if occursin("Finished successfully: optimality conditions satisfied", line)
            converged = true
        end
        
        # Extract x0 (initial design vector)
        if startswith(line, "x0: [")
            try
                # Remove "x0: " prefix and parse the array
                array_str = line[5:end]
                x0 = eval(Meta.parse(array_str))
            catch e
                @warn "Failed to parse x0: $e"
            end
        end
        
        # Extract f0 (initial objective value)
        if startswith(line, "f0: ")
            try
                f0 = parse(Float64, line[5:end])
            catch e
                @warn "Failed to parse f0: $e"
            end
        end
        
        # Extract xopt (optimal design vector)
        if startswith(line, "xopt: [")
            try
                # Remove "xopt: " prefix and parse the array
                array_str = line[7:end]
                xopt = eval(Meta.parse(array_str))
            catch e
                @warn "Failed to parse xopt: $e"
            end
        end
        
        # Extract fopt (optimal objective value)
        if startswith(line, "fopt: ")
            try
                fopt = parse(Float64, line[7:end])
            catch e
                @warn "Failed to parse fopt: $e"
            end
        end
    end
    
    return (
        success = success,
        converged = converged,
        x0 = x0,
        f0 = f0,
        xopt = xopt,
        fopt = fopt,
        perturbation_num = perturbation_num,
        start_time = start_time
    )
end

"""
Parse multiple perturbation log files from a directory.
Returns a vector of results from parse_perturbation_log.
"""
function parse_all_perturbation_logs(directory::String)
    log_files = filter(f -> endswith(f, ".log") && occursin("perturbation", f), readdir(directory))
    results = []
    
    for log_file in log_files
        try
            result = parse_perturbation_log(joinpath(directory, log_file))
            push!(results, result)
            println("Parsed $(log_file): Success = $(result.success), fopt = $(result.fopt)")
        catch e
            @warn "Failed to parse $log_file: $e"
        end
    end
    
    # Sort by perturbation number
    sort!(results, by = r -> r.perturbation_num)
    return results
end


