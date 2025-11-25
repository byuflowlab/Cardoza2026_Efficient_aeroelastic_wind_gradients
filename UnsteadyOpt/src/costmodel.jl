#=
Cost model from WISDEM documentation/
https://wisdem.readthedocs.io/en/master/wisdem/nrelcsm/theory.html

7/22/25
Adam Cardoza
=#

"""
    cost_of_energy(aep, bos, tcc, llc, lrc, o_m; fcr=0.1158, tax_rate=0.4, turbine_number=1)

Calculate the cost of energy (CoE) for a wind turbine.

*Arguments*
- `aep`: Annual energy production in Wh. Todo: Should this be in kWh? -> It doesn't look like it. 
- `bos`: Balance of station costs
- `tcc`: Total capital costs (Cost of the turbine, tower, etc.)
- `llc`: Land lease costs
- `lrc`: Levelized replacement cost
- `o_m`: Operating and maintenance costs

*Keyword Arguments*
- `fcr`: Fixed charge rate, default is 0.1158.
- `tax_rate`: Tax rate, default is 0.4.
- `turbine_number`: Number of turbines, default is 1.

*Returns*
- The cost of energy in USD/kWh... I think

"""
function cost_of_energy(aep, bos, tcc, llc, lrc, o_m; fcr=0.1158, tax_rate=0.4, turbine_number=1)
    warrantyPremium = 0.15*tcc*turbine_number/1.1

    icc = bos + tcc*turbine_number + warrantyPremium #Initial capital costs

    coe = (fcr * (icc))/aep + ((llc + lrc + o_m)*(1 - tax_rate))/aep #WISDEM code equation. 
    return coe
end

function cost_of_energy(aep, bos, tcc, opex; fcr=0.1158, tax_rate=0.4, turbine_number=1)
    warrantyPremium = 0.15*tcc*turbine_number/1.1

    icc = bos + tcc*turbine_number + warrantyPremium #Initial capital costs

    coe = (fcr * (icc))/aep + ((opex)*(1 - tax_rate))/aep #WISDEM code equation. 
    return coe
end

function levelized_cost_of_energy(aep, bos, tcc, llc, lrc, o_m; construction_time = 1, project_lifetime=20., discount_rate=0.07, turbine_number=1)

    warrantyPremium = 0.15*tcc*turbine_number/1.1
    icc = bos + tcc*turbine_number + warrantyPremium

    opex = llc + lrc + o_m

    amortFactor = (1 + ((1 + discount_rate)^ construction_time - 1)) * (discount_rate / (1 - (1 + discount_rate)^(-project_lifetime)))
    lcoe = (icc * amortFactor + opex)/aep

    return lcoe
end

function estimate_blade_mass(D_rotor, class::Int, carbon_spar::Bool; b=nothing)
    k_m = 0.5 #2015 update. 

    if isnothing(b)
        if class == 1 && carbon_spar
            b = 2.47
        elseif class >= 2 && carbon_spar
            b = 2.44
        elseif class == 1 && !carbon_spar
            b = 2.54
        elseif class >= 2 && !carbon_spar
            b = 2.50
        end
    end

    m_blade = k_m*(D_rotor/2)^b 

    return m_blade
end

function blade_cost(m_blade)
    k_c = 14.6

    c_blade = k_c *m_blade #blade cost
    return c_blade
end

function hub_mass(m_blade)
    k_m = 2.3
    b = 1320 #2015 update
    m_hub = k_m * m_blade + b
    return m_hub
end

function hub_cost(m_blade; m_hub=nothing)
    k_c = 3.9
    if isnothing(m_hub)
        m_hub = hub_mass(m_blade)
    end

    c_hub = k_c * m_hub
    return c_hub
end

function bearing_mass(n_blade, m_blade)
    k_m = 0.1295 #2015 update
    b_1 = 491.31
    m_bearing = n_blade * k_m * m_blade + b_1
    return m_bearing
end

function pitch_mass(m_bearing)
    b_2 = 555
    h = 0.328
    m_pitch = m_bearing*(1+h) + b_2
    return m_pitch
end

function pitch_cost(m_pitch)
    k_c = 22.1

    c_pitch = k_c * m_pitch
    return c_pitch
end

function pitch_cost(n_blade::Int, m_blade)
    m_bearing = bearing_mass(n_blade, m_blade)
    m_pitch = pitch_mass(m_bearing)
    c_pitch = pitch_cost(m_pitch)
    return c_pitch
end

function spinner_mass(D_rotor)
    k_m = 15.5 #2015 cost
    b = -980
    m_spin = k_m*D_rotor + b
    return m_spin
end

function spinner_cost(D_rotor; m_spin=nothing)
    if isnothing(m_spin)
        m_spin = spinner_mass(D_rotor)
    end

    k_c = 11.1
    c_spin = k_c * m_spin
    return c_spin
end

function low_speed_shaft_mass(m_blade, P_turbine)
    #P_turbine is the machine rating
    k_m = 13 #2015 update
    b1 = 0.65
    b2 = 775

    # k_m = 20 #On a graph in the theory doc. 
    # b1 = 0.65
    # b2 = -765

    m_lss = (k_m * ((m_blade * P_turbine*1e-3)^b1)) + b2
    return m_lss
end

function low_speed_shaft_cost(m_lss)
    k_c = 11.9
    c_lss = k_c * m_lss
    return c_lss
end

function low_speed_shaft_cost(m_blade, P_turbine)
    #P_turbine is the machine rating
    m_lss = low_speed_shaft_mass(m_blade, P_turbine)
    c_lss = low_speed_shaft_cost(m_lss)
    return c_lss
end

function low_speed_shaft_mass(D_rotor, m_rotor, Q_rotor) #csm_original
    lenShaft = 0.03 * D_rotor
    mmtArm = lenShaft/ 5
    bendLoad = 1.25 * 9.81 * m_rotor
    bendMom = bendLoad * mmtArm
    hFact = 0.1
    hollow = 1/(1-(hFact^4))
    # outDiam = ((32/pi)*hollow*3.25*((Q_rotor*3/371000000.0)^2) + (bendMom/71070000)^2)^(1/3)
    outDiam = ((32.0 / pi) * hollow * 3.25 * ((Q_rotor * 3.0 / 371000000.0)^2 + (bendMom / 71070000)^2)^0.5)^(1.0 / 3.0)
    inDiam = outDiam*hFact

    m_lss = 1.25*(pi/4)*(outDiam^2 - inDiam^2)*lenShaft*7860 
    return m_lss
end

function main_bearings_mass(D_rotor, n_bearing::Int)
    k_m = 0.0001 #2015 update
    b = 3.5

    m_bearing = n_bearing * k_m * D_rotor^b
    return m_bearing
end

function main_bearings_cost(m_bearing)
    k_c = 4.5 #2015 update
    c_bearing = k_c * m_bearing
    return c_bearing
end

function main_bearings_cost(D_rotor, n_bearing::Int)
    m_bearing = main_bearings_mass(D_rotor, n_bearing)
    c_bearing = main_bearings_cost(m_bearing)
    return c_bearing
end

function estimate_rotor_torque(P_turbine, D_rotor, V_tip; eta=0.9)
    #Q_rotor is the rotor torque Q ~ P_turbine * D_rotor/(2*eta*V_tip) where eta is the drivetrain efficiency. 
    Q_rotor = P_turbine * D_rotor / (2 * eta * V_tip)
    return Q_rotor
end

function gearbox_mass(Q_rotor; torque_density=200, update2024::Bool=false)

    if update2024
        m_gearbox = Q_rotor*1e3/torque_density #2024 update
    else
        k_m = 113
        b = 0.71

        m_gearbox = k_m * Q_rotor^b
    end
    
    return m_gearbox
end

function gearbox_cost(Q_rotor; torque_density=200, gearbox_torque_cost=50, m_gearbox=nothing, update2024::Bool=false)
    if isnothing(m_gearbox)
        m_gearbox = gearbox_mass(Q_rotor)
    end

    if update2024
        c_gearbox = torque_density * gearbox_torque_cost * m_gearbox *1e-3 #2024 update
    else
        k_c = 12.9
        c_gearbox = k_c * m_gearbox #Theory doc
    end
    return c_gearbox
end

function brake_mass(Q_rotor)
    # k_m = 1.22 #Theory doc
    k_m = 0.00122 #2015 update (code)
    m_brake = k_m * Q_rotor
    return m_brake
end

function brake_cost(Q_rotor; m_brake=nothing) 
    if isnothing(m_brake)
        m_brake = brake_mass(Q_rotor)
    end
    k_c = 3.6254

    c_brake = k_c * m_brake #2020 update
    return c_brake
end

function high_speed_shaft_mass(P_turbine)
    #P_turbine is the machine rating
    # k_m = 198.94 #Theory doc
    k_m = 0.19894 #2015 update
    m_hss = k_m * P_turbine
    return m_hss
end

function high_speed_shaft_cost(P_turbine; m_hss=nothing)
    #P_turbine is the machine rating

    if isnothing(m_hss)
        m_hss = high_speed_shaft_mass(P_turbine)
    end

    k_c = 6.8 #2015 update
    
    c_hss = k_c * m_hss
    return c_hss
end

function generator_mass(P_turbine)
    #P_turbine is the machine rating

    k_m = 2.3
    b = 3400 #2015 update

    m_gen = k_m * P_turbine + b
    return m_gen
end

function generator_cost(P_turbine; m_gen=nothing)
    #P_turbine is the machine rating
    if isnothing(m_gen)
        m_gen = generator_mass(P_turbine)
    end
    
    k_c = 12.4 #2015 update
    c_gen = k_c * m_gen
    return c_gen
end

function yaw_system_mass(D_rotor)
    #D_rotor is the rotor diameter

    k_m = 0.0009
    b = 3.314
    # m_yaw = k_m * D_rotor^b #Theory doc
    m_yaw = 1.5 * k_m * D_rotor^b #2015 update (code)
    return m_yaw
end

function yaw_system_cost(D_rotor; m_yaw=nothing)
    if isnothing(m_yaw)
        m_yaw = yaw_system_mass(D_rotor)
    end

    k_c = 8.3 #Not updated in 2015
    c_yaw = k_c * m_yaw
    return c_yaw
end

function hydraulic_cooling_mass(P_turbine)
    #P_turbine is the machine rating
    # k_m = 80 #Theory doc
    k_m = 0.08 #2015 update (code)
    m_hvac = k_m * P_turbine
    return m_hvac
end

function hydraulic_cooling_cost(P_turbine; m_hvac=nothing)
    #P_turbine is the machine rating
    if isnothing(m_hvac)
        m_hvac = hydraulic_cooling_mass(P_turbine)
    end

    k_c = 124 #2015 update
    c_hvac = k_c * m_hvac
    return c_hvac
end

function converter_cost(m_converter)
    k_c = 18.8 #2015 update 
    c_converter = k_c * m_converter
    return c_converter
end

function transformer_mass(P_turbine)
    #P_turbine is the machine rating
    k_m = 1.915 #2015 update
    b = 1910

    m_transformer = k_m * P_turbine + b
    return m_transformer
end

function transformer_cost(P_turbine; m_transformer=nothing)
    if isnothing(m_transformer)
        m_transformer = transformer_mass(P_turbine)
    end
    k_c = 18.8 #2015 update 
    c_transformer = k_c * m_transformer
    return c_transformer
end

function connect_cost(P_turbine) #Called electrial connection cost
    # return 41850*P_turbine #Theory doc
    return 41.85 * P_turbine #2015 update (code)
end

function control_cost(P_turbine)
    # return 21150*P_turbine #Theory doc
    return 21.15 * P_turbine #2015 update (code)
end

function platform_mass(m_bedplate; crane::Bool=true)
    k_m = 0.125

    m_platform = k_m * m_bedplate
    if crane #2015 update
        m_platform += crane_mass()
    end

    return m_platform
end

function platform_cost(m_bedplate; crane::Bool=true, m_platform=nothing)
    if isnothing(m_platform)
        m_platform = platform_mass(m_bedplate; crane)
    end

    k_c = 17.1 #2015 update
    c_platform = k_c * m_platform

    if crane #2015 update
        c_platform += crane_cost()
    end
    return c_platform
end

function crane_mass()
    m_crane = 3000 #mass of the crane in kg
    return m_crane
end

function crane_cost()
    c_crane = 12000
    return c_crane
end

function bedplate_mass(D_rotor)
    b = 2.2 #2015 update

    m_bedplate = D_rotor^b
    return m_bedplate
end

function bedplate_cost(D_rotor; m_bedplate=nothing)
    if isnothing(m_bedplate)
        m_bedplate = bedplate_mass(D_rotor)
    end

    k_c = 2.9
    c_bedplate = k_c * m_bedplate
    return c_bedplate
end

function nacelle_cover_mass(P_turbine)
    k_m = 1.2817 #2015 update
    b = 428.19

    m_cover = k_m * P_turbine + b
    return m_cover
end


function nacelle_cover_cost(P_turbine; m_cover=nothing)
    if isnothing(m_cover)
        m_cover = nacelle_cover_mass(P_turbine)
    end
    k_c = 5.7 #2015 update
    c_cover = k_c * m_cover
    return c_cover
end

function tower_mass(L_tower)
    k_m = 19.828
    b = 2.0282

    m_tower = k_m * L_tower^b
    return m_tower
end

function tower_cost(L_tower)
    k_m = 19.828
    b = 2.0282
    k_c = 2.9 #2015 update

    m_tower = k_m * L_tower^b
    c_tower = k_c * m_tower
    return c_tower
end

function hub_system_mass(m_hub, m_pitch, m_spin)
    m_hubsys = m_hub + m_pitch + m_spin
    return m_hubsys
end

function hub_system_cost(c_hub, c_pitch, c_spin; kt_hs=0.0, kp_hs=0.0, ko_hs=0.0, ka_hs=0.0)

    #default kt = kp = ko = ka = 0.0
    #kt is a transport multiplier
    #kp is a profit multiplier
    #ko is an overhead cost multiplier
    #ka is an assembly cost multiplier
    c_hubsys = (1 + kt_hs + kp_hs)*(1 + ko_hs + ka_hs)*(c_hub + c_pitch + c_spin)
    return c_hubsys
end

function rotor_system_mass(n_blade::Int, m_blade, m_hubsys)
    m_rotor = n_blade*m_blade + m_hubsys
    return m_rotor
end

function rotor_system_cost(n_blade::Int, c_blade, c_hubsys)
    c_rotor = n_blade*c_blade + c_hubsys
    return c_rotor
end

function nacelle_mass(m_lss, m_bearing, m_gearbox, m_hss, m_gen, m_bedplate, m_yaw, m_hvac, m_transformer, m_platform, m_cover)
    return m_lss + m_bearing + m_gearbox + m_hss + m_gen + m_bedplate + m_yaw + m_hvac + m_transformer + m_platform + m_cover
end

function nacelle_cost(c_lss, c_bearing, c_gearbox, c_hss, c_brake, c_gen, c_bedplate, c_yaw, c_converter, c_hvac, c_transformer, c_connect, c_control, c_platform, c_cover; kt_na=0.0, kp_na=0.0, k0_na=0.0, ka_na=0.0)
    c_parts = c_lss + c_bearing + c_gearbox + c_hss + c_brake + c_gen + c_bedplate + c_yaw + c_converter + c_hvac + c_transformer + c_connect + c_control + c_platform + c_cover

    c_nacelle = (1 + kt_na + kp_na)*(1 + k0_na + ka_na)*c_parts
    return c_nacelle
end

function tower_system_cost(c_tower; kt_t=0.0, kp_t=0.0, k0_t=0.0, ka_t=0.0)
    c_towersys = (1 + kt_t + kp_t)*(1 + k0_t + ka_t)*c_tower
    return c_towersys
end

function turbine_mass(m_rotor, m_nacelle, m_tower)
    m_turbine = m_rotor + m_nacelle + m_tower
    return m_turbine
end

function turbine_cost(c_rotor, c_nacelle, c_tower; kt_tu=0.0, kp_tu=0.0, k0_tu=0.0, ka_tu=0.0)
    c_turbine = (1 + kt_tu + kp_tu)*(1 + k0_tu + ka_tu)*(c_rotor + c_nacelle + c_tower)
    return c_turbine
end

function calc_turbine_cost(P_turbine, D_rotor, Q_rotor, n_blade, m_blade, n_bearing, L_tower; crane::Bool=true, m_converter=0.0,
    kt_tu=0.0, kp_tu=0.0, k0_tu=0.0, ka_tu=0.0, kt_t=0.0, kp_t=0.0, k0_t=0.0, ka_t=0.0, kt_na=0.0, kp_na=0.0, k0_na=0.0, ka_na=0.0, kt_hs=0.0, kp_hs=0.0, ko_hs=0.0, ka_hs=0.0)

    c_blade = blade_cost(m_blade)
    c_hub = hub_cost(m_blade)
    c_pitch = pitch_cost(n_blade, m_blade)
    c_spin = spinner_cost(D_rotor)
    c_hubsys = hub_system_cost(c_hub, c_pitch, c_spin; kt_hs, kp_hs, ko_hs, ka_hs)
    c_rotor = rotor_system_cost(n_blade, c_blade, c_hubsys)


    c_lss = low_speed_shaft_cost(m_blade, P_turbine)
    c_bearing = main_bearings_cost(D_rotor, n_bearing)
    c_gearbox = gearbox_cost(Q_rotor) 
    c_hss = high_speed_shaft_cost(P_turbine)
    c_brake = brake_cost(Q_rotor)
    c_gen = generator_cost(P_turbine)
    m_bedplate = bedplate_mass(D_rotor) 
    c_bedplate = bedplate_cost(D_rotor) 
    c_yaw = yaw_system_cost(D_rotor)
    c_converter = converter_cost(m_converter)
    c_hvac = hydraulic_cooling_cost(P_turbine)
    c_transformer = transformer_cost(P_turbine)
    c_connect = connect_cost(P_turbine)
    c_control = control_cost(P_turbine)
    c_platform = platform_cost(m_bedplate; crane)
    c_cover = nacelle_cover_cost(P_turbine)
    c_nacelle = nacelle_cost(c_lss, c_bearing, c_gearbox, c_hss, c_brake, c_gen, c_bedplate, c_yaw, c_converter, c_hvac, c_transformer, c_connect, c_control, c_platform, c_cover; kt_na, kp_na, k0_na, ka_na)


    c_tower = tower_cost(L_tower)
    c_towersys = tower_system_cost(c_tower; kt_t, kp_t, k0_t, ka_t)

    return turbine_cost(c_rotor, c_nacelle, c_towersys; kt_tu, kp_tu, k0_tu, ka_tu)
end


function balance_of_station(machine_rating, rotor_diameter, hub_height; bos_multiplier=1.0, turbine_number=1, multiplier=1.0, ref_year=2002, ref_month=9, end_year=2009, end_month=12)

    lPrmtsCostCoeff1 = 9.94e-04
    lPrmtsCostCoeff2 = 20.31
    oPrmtsCostFactor = 37.0  # $/kW (2003)
    scourCostFactor = 55.0  # $/kW (2003)
    ptstgCostFactor = 20.0  # $/kW (2003)
    ossElCostFactor = 260.0  # $/kW (2003) shallow
    ostElCostFactor = 290.0  # $/kW (2003) transitional
    ostSTransFactor = 25.0  # $/kW (2003)
    ostTTransFactor = 77.0  # $/kW (2003)
    osInstallFactor = 100.0  # $/kW (2003) shallow & trans
    suppInstallFactor = 330.0  # $/kW (2003) trans additional
    paiCost = 60000.0  # per turbine

    suretyBRate = 0.03  # 3% of ICC
    suretyBond = 0.0
    
    fcCoeff = 303.23
    fcExp = 0.4037

    SweptArea = ((rotor_diameter/2)^2)*pi
    foundation_cost = fcCoeff * (hub_height * SweptArea)^fcExp
    fndnCostEscalator = IPPI_FND(ref_year, ref_month, end_year, end_month) 
    foundation_cost *= fndnCostEscalator

    # self.d_foundation_d_diameter = (fndnCostEscalator * fcCoeff* fcExp
    #     * ((hub_height * (2.0 * 0.5 * (rotor_diameter * 0.5) * pi)) ^ (fcExp - 1))
    #     * hub_height)

    # self.d_foundation_d_hheight = (fndnCostEscalator * fcCoeff * fcExp * ((self.hub_height * SweptArea) ** (fcExp - 1)) * SweptArea)


    # cost calculations
    tpC1 = 0.00001581
    tpC2 = -0.0375
    tpInt = 54.7
    tFact = (tpC1 * (machine_rating^2)) + (tpC2*machine_rating) + tpInt

    roadsCivil_costs = 0.0
    portStaging_costs = 0.0
    pai_costs = 0.0
    scour_costs = 0.0


    engPermits_costs = (lPrmtsCostCoeff1 * machine_rating * machine_rating) + (lPrmtsCostCoeff2 * machine_rating)

    ref_month = 3 #Todo: Who knows why they keep changing this.... 

    engPermits_costs *= IPPI_LPM(ref_year, ref_month, end_year, end_month)

    # self.d_development_d_rating = ppi.compute("IPPI_LPM") * (
    #     2.0 * lPrmtsCostCoeff1 * self.machine_rating + lPrmtsCostCoeff2)

    ref_month = 9 #Todo: Who knows why they keep changing this.... 

    elC1 = 3.49e-06
    elC2 = -0.0221
    elInt = 109.7
    eFact = elC1 * machine_rating * machine_rating + elC2 * machine_rating + elInt
    electrical_costs = machine_rating * eFact * IPPI_LEL(ref_year, ref_month, end_year, end_month)

    # self.d_electrical_d_rating = ppi.compute("IPPI_LEL") * (
    #     3.0 * elC1 * self.machine_rating**2.0 + 2.0 * elC2 * self.machine_rating + elInt)

    rcC1 = 2.17e-06
    rcC2 = -0.0145
    rcInt = 69.54
    rFact = rcC1 * machine_rating * machine_rating + rcC2 * machine_rating + rcInt
    roadsCivil_costs = machine_rating * rFact * IPPI_RDC(ref_year, ref_month, end_year, end_month)
    # self.d_preparation_d_rating = ppi.compute("IPPI_RDC") * (
    #     3.0 * rcC1 * self.machine_rating**2.0 + 2.0 * rcC2 * self.machine_rating + rcInt)

    iCoeff = 1.965
    iExp = 1.1736
    installation_costs = iCoeff * ((hub_height * rotor_diameter)^iExp) * IPPI_LAI(ref_year, ref_month, end_year, end_month)
    # self.d_assembly_d_diameter = (iCoeff
    #     * ((self.hub_height * self.rotor_diameter) ** (iExp - 1))
    #     * self.hub_height
    #     * ppi.compute("IPPI_LAI"))
    # self.d_assembly_d_hheight = (iCoeff
    #     * ((self.hub_height * self.rotor_diameter) ** (iExp - 1))
    #     * self.rotor_diameter
    #     * ppi.compute("IPPI_LAI"))

    transportation_costs = machine_rating * tFact * IPPI_TPT(ref_year, ref_month, end_year, end_month)

    # self.d_transport_d_rating = ppi.compute("IPPI_TPT") * (
        # tpC1 * 3.0 * self.machine_rating**2.0 + tpC2 * 2.0 * self.machine_rating + tpInt)


    bos_costs = foundation_cost + transportation_costs + roadsCivil_costs + portStaging_costs
        + installation_costs + electrical_costs + engPermits_costs + pai_costs + scour_costs




    # self.d_other_d_tcc = 0.0


    suretyBond = 0.0


    bos_costs = turbine_number * (bos_costs + suretyBond)
    bos_costs *= bos_multiplier  

    # self.bos_breakdown_development_costs = engPermits_costs * self.turbine_number
    # self.bos_breakdown_preparation_and_staging_costs = (roadsCivil_costs + portStaging_costs) * self.turbine_number
    # self.bos_breakdown_transportation_costs = transportation_costs * self.turbine_number
    # self.bos_breakdown_foundation_and_substructure_costs = foundation_cost * self.turbine_number
    # self.bos_breakdown_electrical_costs = electrical_costs * self.turbine_number
    # self.bos_breakdown_assembly_and_installation_costs = installation_costs * self.turbine_number
    # self.bos_breakdown_soft_costs = 0.0
    # self.bos_breakdown_other_costs = (pai_costs + scour_costs + suretyBond) * self.turbine_number

    return bos_costs
end

function operating_maintenance_costs(aep; ref_year=2002, ref_month=9, end_year=2009, end_month=12)

    landCostFactor = 0.0070  # $/kwH

    cost = aep * landCostFactor
    costEscalator = IPPI_LOM(ref_year, ref_month, end_year, end_month)

    return cost * costEscalator
end

function levelized_replacement_costs(machine_rating; ref_year=2002, ref_month=9, end_year=2009, end_month=12, turbine_number=1)
    lrcCF = 10.70  # land based
    costlrcEscFactor = IPPI_LLR(ref_year, ref_month, end_year, end_month)

    return machine_rating * lrcCF * costlrcEscFactor * turbine_number
end

function land_lease_costs(aep; ref_year=2002, ref_month=9, end_year=2009, end_month=12)
    leaseCF = 0.00108  # land based
    costlandEscFactor = IPPI_LSE(ref_year, ref_month, end_year, end_month)

    return aep * leaseCF * costlandEscFactor
end