"""
    read_precomp_sections_file(file, loc_w)

Reads composite sections from a PreComp input file
"""
function read_precomp_sections_file(file, loc_w)

    open(file) do f
        readline(f)
        readline(f)
        readline(f)

        # number of sectors
        nu = parse(Int, split(chomp(readline(f)))[1])

        readline(f)
        readline(f)

        # read normalized chord locations
        loc_u = [parse(Float64, x) for x in split(chomp(readline(f)))]

        nlam_u, nply_u, tply_u, theta_u, imat_u = read_precomp_sector_from_file(f, nu)

        readline(f)
        readline(f)
        readline(f)

        # number of sectors
        nl = parse(Int, split(chomp(readline(f)))[1])

        readline(f)
        readline(f)

        loc_l = [parse(Float64, x) for x in split(chomp(readline(f)))]

        nlam_l, nply_l, tply_l, theta_l, imat_l = read_precomp_sector_from_file(f, nl)

        readline(f)
        readline(f)
        readline(f)
        readline(f)

        nw = length(loc_w)

        nlam_w, nply_w, tply_w, theta_w, imat_w = read_precomp_sector_from_file(f, nw)

        return loc_u, nlam_u, nply_u, tply_u, theta_u, imat_u,
            loc_l, nlam_l, nply_l, tply_l, theta_l, imat_l,
            loc_w, nlam_w, nply_w, tply_w, theta_w, imat_w
    end
end

"""
read_precomp_sector_from_file(f, n)

Reads PreComp sector from a precomp file.
"""
function read_precomp_sector_from_file(f, n)

    nlam = Vector{Int}(undef, n)
    nply = Vector{Vector{Int}}(undef, n)
    tply = Vector{Vector{Float64}}(undef, n)
    theta = Vector{Vector{Float64}}(undef, n)
    imat = Vector{Vector{Int}}(undef, n)

    for i = 1:n
        readline(f)
        readline(f)

        line = chomp(readline(f))

        nlam[i] = parse(Int, split(line)[2])

        readline(f)
        readline(f)
        readline(f)
        readline(f)

        nply[i] = zeros(Int, nlam[i])
        tply[i] = zeros(Float64, nlam[i])
        theta[i] = zeros(Float64, nlam[i])
        imat[i] = zeros(Int, nlam[i])

        for j = 1:nlam[i]
            array = split(chomp(readline(f)))
            nply[i][j] = parse(Int, array[2])
            tply[i][j] = parse(Float64, array[3])
            theta[i][j] = parse(Float64, array[4])*pi/180
            imat[i][j] = parse(Int, array[5])
        end
    end

    return nlam, flatten(nply), flatten(tply), flatten(theta), flatten(imat)
end

flatten(x) = collect(Iterators.flatten(x))

"""
    read_precomp_profile_file(file)

Reads a precomp profile
"""
function read_precomp_profile_file(file)
    open(file, "r") do f

        # read header
        naf = parse(Int, split(chomp(readline(f)))[1])
        readline(f)
        readline(f)
        readline(f)

        # read airfoil data
        x = zeros(naf)
        y = zeros(naf)
        for i = 1:naf
            data = split(chomp(readline(f)))
            x[i] = parse(Float64, data[1])
            y[i] = parse(Float64, data[2])
        end

        return x, y
    end
end

"""
    read_precomp_materials_file(file)

reads material properties from PreComp input file
"""
function read_precomp_materials_file(file)

    open(file) do f

        # skip through header
        for i in 1:3
            readline(f)
        end

        mat_id = Int[]
        E1 = Float64[]
        E2 = Float64[]
        G12 = Float64[]
        nu12 = Float64[]
        rho = Float64[]
        name = String[]
        for line in eachline(f)
            array = split(line, limit=7)
            push!(mat_id, parse(Int, array[1]))
            push!(E1, parse(Float64, array[2]))
            push!(E2, parse(Float64, array[3]))
            push!(G12, parse(Float64, array[4]))
            push!(nu12, parse(Float64, array[5]))
            push!(rho, parse(Float64, array[6]))
            push!(name, strip(chomp(array[7])))
        end

        isort = sortperm(mat_id)
        E1 = E1[isort]
        E2 = E2[isort]
        G12 = G12[isort]
        nu12 = nu12[isort]
        rho = rho[isort]
        name = name[isort]

        return E1, E2, G12, nu12, rho, name
    end
end

"""
    get_precomp_descriptions(path)
Read in the precomp files and return the descriptions. Skip the first and last point because they are close to the hub and the tip. 
"""
function get_precomp_descriptions(path, yamlpath)
    file = YAML.load_file(joinpath(yamlpath, "distro.yaml"))

    rvec = file["radius"][2:end-1]
    cvec = file["chord"][2:end-1]
    twistvec = file["twist"][2:end-1].*(pi/180)
    le_loc = file["le_str"][2:end-1]

    distro = (; rvec, cvec, twistvec, le_loc)

    ### Read in the materials and convert them. 
    e1, e2, g12, nu12, rho, names = read_precomp_materials_file(joinpath(path,"materials.inp"))

    gelcoat = GXBeamCS.Material(e1[1], e2[1], e2[1], g12[1], g12[1], g12[1], nu12[1], nu12[1], nu12[1], rho[1], e1[1]/100, e1[1]/100, e2[1]/100, e2[1]/100, e2[1]/100, e2[1]/100, g12[1]/100, g12[1]/100, g12[1]/100)

    ELT = GXBeamCS.Material(e1[2], e2[2], e2[2], g12[2], g12[2], g12[2], nu12[2], nu12[2], nu12[2], rho[2], e1[2]/100, e1[2]/100, e2[2]/100, e2[2]/100, e2[2]/100, e2[2]/100, g12[2]/100, g12[2]/100, g12[2]/100)

    Triax = GXBeamCS.Material(e1[3], e2[3], e2[3], g12[3], g12[3], g12[3], nu12[3], nu12[3], nu12[3], rho[3], e1[3]/100, e1[3]/100, e2[3]/100, e2[3]/100, e2[3]/100, e2[3]/100, g12[3]/100, g12[3]/100, g12[3]/100)

    Saertex = GXBeamCS.Material(e1[4], e2[4], e2[4], g12[4], g12[4], g12[4], nu12[4], nu12[4], nu12[4], rho[4], e1[4]/100, e1[4]/100, e2[4]/100, e2[4]/100, e2[4]/100, e2[4]/100, g12[4]/100, g12[4]/100, g12[4]/100)

    foam = GXBeamCS.Material(e1[5], e2[5], e2[5], g12[5], g12[5], g12[5], nu12[5], nu12[5], nu12[5], rho[5], e1[5]/100, e1[5]/100, e2[5]/100, e2[5]/100, e2[5]/100, e2[5]/100, g12[5]/100, g12[5]/100, g12[5]/100)

    carbon = GXBeamCS.Material(e1[6], e2[6], e2[6], g12[6], g12[6], g12[6], nu12[6], nu12[6], nu12[6], rho[6], e1[6]/100, e1[6]/100, e2[6]/100, e2[6]/100, e2[6]/100, e2[6]/100, g12[6]/100, g12[6]/100, g12[6]/100)

    materials = [gelcoat, ELT, Triax, Saertex, foam, carbon]




    ### Running through each station
    w1vec = file["web1"]
    w2vec = file["web2"]
    w3vec = file["web3"]

    stations = ()

    for k in eachindex(rvec)
        fn = k+1
        ### Accessing airfoil data 
        x, y = read_precomp_profile_file(joinpath(path, "shape_"*string(fn)*".inp"))
        # Rearrange the order.
        xmaxidx = argmax(x)
        xtop = reverse(x[1:xmaxidx])
        ytop = reverse(y[1:xmaxidx])
        xbot = reverse(x[xmaxidx:end])
        ybot = reverse(y[xmaxidx:end])
        xaf = vcat(xtop, xbot)
        yaf = vcat(ytop, ybot)


        ### Accessing the layups for this section
        # Need the web locations apparently. 
        w1 = w1vec[fn] >= 0 ? w1vec[fn] : missing
        w2 = w2vec[fn] >= 0 ? w2vec[fn] : missing
        w3 = w3vec[fn] >= 0 ? w3vec[fn] : missing
        webloc = collect(skipmissing([w1, w2, w3]))

        locU, n_laminaU, n_pliesU, tU, thetaU, mat_idxU,
        locL, n_laminaL, n_pliesL, tL, thetaL, mat_idxL,
        locW, n_laminaW, n_pliesW, tW, thetaW, mat_idxW =
        read_precomp_sections_file(joinpath(path, "layup_"*string(fn)*".inp"), webloc)

        ## Get the layup and convert it to GXBeamCS layers
        #### Layups
        nl = length(n_laminaU)
        total_idxs = zeros(Int, nl+1)
        total_idxs[2:end] = [sum(n_laminaU[1:i]) for i in 1:length(n_laminaU)]
        segments = Vector{Vector{GXBeamCS.Layer{Float64}}}(undef, nl)
        layer_idxs = Vector{Vector{Int}}(undef, nl)
        for i in 1:nl
            idxs = total_idxs[i]+1:total_idxs[i+1]
            matidx = mat_idxU[idxs]  # material index
            t = tU[idxs]
            theta = thetaU[idxs]
            segments[i] = GXBeamCS.Layer.(materials[matidx], t, theta)
            layer_idxs[i] = matidx
        end

        xbreak = locU

        #### Webs
        web_idxs = [sum(n_laminaW[1:i]) for i in 1:length(n_laminaW)]
        nwebs = 0
        nwebs += ismissing(w1) ? 0 : 1
        nwebs += ismissing(w2) ? 0 : 1
        nwebs += ismissing(w3) ? 0 : 1

        webs_layer_idxs = Vector{Vector{Int}}(undef, nwebs)

        if ismissing(w1)
            web1 = missing
        else
            widxs = 1:web_idxs[1]
            matidx = mat_idxW[widxs]  # material index
            t = tW[widxs]
            theta = thetaW[widxs]
            web1 = GXBeamCS.Layer.(materials[matidx], t, theta)
            webs_layer_idxs[1] = matidx
        end


        if ismissing(w2)
            web2 = missing
        else
            widxs = web_idxs[1]:web_idxs[2]
            matidx = mat_idxW[widxs]  # material index
            t = tW[widxs]
            theta = thetaW[widxs]
            matidx = vcat(mat_idxW[1], matidx[3:4])
            t = vcat(sum(t[1:2]), t[3:4])
            theta = vcat(sum(theta[1:2]), theta[3:4])
            web2 = GXBeamCS.Layer.(materials[matidx], t, theta)
            webs_layer_idxs[2] = matidx
        end

        if ismissing(w3)
            web3 = missing
        else
            web3 = missing
            webloc = webloc[1:2] 
        end

         
        webs = collect(skipmissing([web1, web2, web3]))
        
        stations = (stations..., (;xaf, yaf, xbreak, webloc, segments, webs, layer_idxs, webs_layer_idxs))
    end

    return distro, materials, stations
end



### Cross section information (of the modified 5-segment formulation)
num_segments5 = ones(Int, 36).*5
Nsegs5::Int64 = sum(Int, num_segments5)
num_elements5 = [632, 512, 512, 512, 512, 512, 644, 658, 658, 658, 658, 658, 528, 750, 580, 748, 580, 576, 576, 746, 574, 746, 744, 574, 742, 574, 736, 578, 578, 392, 526, 356, 356, 356, 356, 312] 
N_elements5::Int64 = sum(Int, num_elements5)
N_bulk_full5::Int64 = 2*N_elements5

function get_scaling_factor(fvec, fn, i, ns)
    if fn<2
        return fvec[i]
    else
        return fvec[sum(Int, ns[1:fn-1]) + i]
    end
end


function scale_segments(stations, fvec, fw1, fw2, materials, num_segs)

    TF = promote_type(typeof(fvec[1]), typeof(stations[1].segments[1][1].t))

    new_stations = ()
    for fn in eachindex(stations)
        
        #### Layups
        ns = length(stations[fn].segments)
        
        segments = Vector{Vector{GXBeamCS.Layer{TF}}}(undef, ns)

        for i in 1:ns #Iterate across the airfoil segments
            layers = stations[fn].segments[i]

            f_thick = get_scaling_factor(fvec, fn, i, num_segs)

            #Scale the layer thicknesses of this section
            t = [l.t for l in layers].*f_thick
            
            materials_i = materials[stations[fn].layer_idxs[i]]
            theta = [(l.theta) for l in layers]

            segments[i] = GXBeamCS.Layer.(materials_i, t, theta)
        end

        #### Webs
        nw = length(stations[fn].webs)
        webs = Vector{Vector{GXBeamCS.Layer{TF}}}(undef, nw)

        for i in 1:nw #Iterate across the web segments
            layers = stations[fn].webs[i]

            if fn > 6
                if i == 1
                    f_thick = fw1[fn-6]
                else
                    f_thick = fw2[fn-6]
                end
            else
                f_thick = 1
            end

            t = [l.t for l in layers].*f_thick

            materials_i = materials[stations[fn].webs_layer_idxs[i]]
            theta = [TF(l.theta) for l in layers]

            webs[i] = GXBeamCS.Layer.(materials_i, t, theta)
        end
        
        new_stations = (new_stations..., (; segments, webs))
    end

    return new_stations
end

function scale_segments(stations, fvec, materials, num_segs)

    TF = promote_type(typeof(fvec[1]), typeof(stations[1].segments[1][1].t))

    new_stations = ()
    for fn in eachindex(stations)
        

        #### Layups
        ns = length(stations[fn].segments)
        
        segments = Vector{Vector{GXBeamCS.Layer{TF}}}(undef, ns)

        for i in 1:ns #Iterate across the airfoil segments
            layers = stations[fn].segments[i]

            f_thick = get_scaling_factor(fvec, fn, i, num_segs)

            #Scale the layer thicknesses of this section
            t = [l.t for l in layers].*f_thick
            
            materials_i = materials[stations[fn].layer_idxs[i]]
            theta = [TF(l.theta) for l in layers]

            segments[i] = GXBeamCS.Layer.(materials_i, t, theta)
        end

        #### Webs
        nw = length(stations[fn].webs)
        webs = Vector{Vector{GXBeamCS.Layer{TF}}}(undef, nw)

        for i in 1:nw #Iterate across the web segments
            layers = stations[fn].webs[i]

            f_thick = TF(1) 

            t = [l.t for l in layers].*f_thick

            materials_i = materials[stations[fn].webs_layer_idxs[i]]
            theta = [TF(l.theta) for l in layers]

            webs[i] = GXBeamCS.Layer.(materials_i, t, theta)
        end
        
        new_stations = (new_stations..., (; segments, webs))
    end

    return new_stations
end


"""
    get_clt()
Create a CLT object for the optimization. 
"""
function get_clt(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs)
    sections = GXBeamCS.get_beam_sections(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs) #TODO: Potentially could directly augment the sections? 
    closed_section = true

    return GXBeamCS.CLT(sections, closed_section)
end

function get_clt_sections_oop(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs)
    sections = GXBeamCS.get_beam_sections_oop(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs) #TODO: Potentially could directly augment the sections? 

    return sections
end



########## Failure constraint functions ############

function get_failure_constraint_indices(i, num_buckling, num_strain, N_buckling)
    if i == 1
        b_idx = 1:num_buckling[i]
        fs_idx = N_buckling+1:N_buckling+num_strain[i]
    else
        b_idx = sum(Int, num_buckling[1:i-1])+1:sum(Int, num_buckling[1:i])
        fs_idx = N_buckling+sum(Int, num_strain[1:i-1])+1:N_buckling+sum(Int, num_strain[1:i])
    end
    return b_idx, fs_idx
end

const num_buckling5 = Int64[20, 20, 20, 20, 20, 20, 20, 25, 20, 25, 20, 20, 16, 20, 16, 20, 16, 24, 24, 28, 24, 28, 28, 24, 28, 28, 36, 32, 32, 20, 32, 20, 20, 20, 20, 36]
const N_buckling5 = sum(Int, num_buckling5)

function check_failure(F, M, clt, strain_ult; top::Bool=true) 
    _, _, strain_p, _ = strains_and_stresses(F, M, clt) 

   if top
        N = size(strain_p, 2)
        strain = view(strain_p, 1, 1:2:N) #Top of each cell
        spar_cap_indices = GXBeamCS.get_section_indices(clt, 3) #Spar cap
        spar_cap_strain = view(strain, spar_cap_indices) #Indexing the strain with the section indices (note not the strain indices which correspond to the top and bottom of each cell). This works because I skip either the top or the bottom of each cell. 
    else
        strain = view(strain_p, 1, 2:2:N) #Bottom of each cell
        spar_cap_indices = GXBeamCS.get_section_indices(clt, 3) #Spar cap
        spar_cap_strain = view(strain, spar_cap_indices)
    end

    ### Buckling
    section = clt.sections[3] #Spar cap 
    b = section.y[1] - section.y[end]

    z, _ = GXBeamCS.zspacing(section.laminate)
    h = z[end] - z[1]
    
    A, B, D = GXBeamCS.laminatestiffnessmatrix(section.laminate, z)

    S = [A B; B D]
    E_axial = det(S)/det(S[2:end, 2:end])/h

    dterm = sqrt(D[1, 1]*D[2, 2]) + D[1, 2] + 2*D[3,3]
    Ncrit = 3.6*((pi/b)^2)*dterm

    eps_crit = -Ncrit/(E_axial*h)

    bm = @. spar_cap_strain*1.35 - eps_crit #buckling margin

    fb = bm./strain_ult


    ### Ultimate strain
    fs = strain./strain_ult #Axial strain
    
    return fb, fs, eps_crit
end






############ Fatigue constraint Function  #################
"""
    get_fatigue_constaint_indces(fat_idxs, j)
Get the indices of the fatigue constraints for a given station. Note: These will need to be shifted by the total number of other constraints. 

**Arguments**
- `fat_idxs`: The radial stations at which we'll calculate fatigue constraints.
- `j`: The index of the station we want to get the indices for.

**Returns**
- `d_idx`: The indices to be shifted of the fatigue constraints for the given station. (They need to be shifted by all the other constraints.)
"""
function get_fatigue_constraint_indices(fat_idxs, j, num_elem)
    if j == 1
        d_idx = 1:num_elem[fat_idxs][1]

    else
        d_idx = sum(Int, num_elem[fat_idxs][1:j-1])+1:sum(Int, num_elem[fat_idxs][1:j])
    end

    return d_idx
end


mean(x) = sum(x)/length(x)

function calculate_damage(g, fat_idxs, forces, moments, clt_list, eps_ult, m, t_elapsed, nu, num_constraints, num_elem; years20 = 20*365*24*60*60, top::Bool=true)
    # This damage is going to be 
    nt = size(forces)[2]
    TF = typeof(forces[1, 1, 1])

    Ttot = years20/t_elapsed
    
    for k in eachindex(fat_idxs)
        num_cells2_k = num_elem[fat_idxs[k]] #The number of strains in the cross section (top and bottom of each cell). 
        axial_strains = zeros(TF, nt-200, num_cells2_k) #Only going to look at the axial strains because they should be the largest (and consistently have been).  
        for i in 201:nt
            fi = forces[k, i, :]
            mi = moments[k, i, :] 

            ### Calculate the damage equivalent strains
            _, _, strains, _ = strains_and_stresses(fi, mi, clt_list[fat_idxs[k]]) 
            N = size(strains, 2) #Number of cells in the cross section
            if top
                axial_strains[i-200, :] = view(strains, 1, 1:2:N) #Top of each cell
            else
                axial_strains[i-200, :] = view(strains, 1, 2:2:N) #Bottom of each cell
            end
        end

        didx = get_fatigue_constraint_indices(fat_idxs, k, num_elem)
        didx = didx .+ num_constraints #shift didx by the number of other constraints.

        for i = 1:num_cells2_k
            Dj, _ , _ = of.damage(axial_strains[:, i]; m=m, Lult=eps_ult)
            g[didx[i]] = log(Dj*Ttot/nu)/200 #Dlife and scaling factor
        end
    end
end

function calculate_damage_at_index(num_elem, fat_idxs, k, forces, moments, clt_list, eps_ult, m; top::Bool=true)
    nt = size(forces)[2]
    TF = typeof(forces[1, 1, 1])
    
    num_cells2_k = num_elem[fat_idxs[k]] #The number of strains in the cross section (top and bottom of each cell). 
    axial_strains = zeros(TF, nt-200, num_cells2_k) #Only going to look at the axial strains because they should be the largest (and consistently have been).  

    for i in 201:nt
        fi = forces[k, i, :]
        mi = moments[k, i, :] 

        ### Calculate the damage equivalent strains
        _, _, strains, _ = strains_and_stresses(fi, mi, clt_list[fat_idxs[k]]) 

        N = size(strains, 2) #Number of cells in the cross section
        if top
            axial_strains[i-200, :] = view(strains, 1, 1:2:N) #Top of each cell
        else
            axial_strains[i-200, :] = view(strains, 1, 2:2:N) #Bottom of each cell
        end
    end

    return [of.damage(axial_strains[:, i]; m=m, Lult=eps_ult)[1] for i in 1:num_cells2_k] #simulation damage values
end

function calculate_damage_oop(fat_idxs, forces, moments, clt_list, eps_ult, m, t_elapsed, num_elem; years20 = 20*365*24*60*60, top::Bool=true)

    Ttot = years20/t_elapsed
    
    damages = [calculate_damage_at_index(num_elem, fat_idxs, k, forces, moments, clt_list, eps_ult, m; top=top) for k in eachindex(fat_idxs)]

    D = mapreduce(identity, vcat, damages).*Ttot

    return D
end

function get_damage(g, fat_idxs, num_elements, num_constraints)
    nD = sum(Int, num_elements[fat_idxs])
    D = zeros(nD)

    for k in eachindex(fat_idxs)
        didx = get_fatigue_constraint_indices(fat_idxs, k, num_elements)
        D[didx] = g[didx.+num_constraints]
    end

    return exp.(D.*200)
end
