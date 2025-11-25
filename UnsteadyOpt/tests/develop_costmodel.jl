#=
Developing the cost model and making sure the functions work as expected.

8/9/25 Adam Cardoza

#Todo: The docs have example masses for the NREL 5MW so I can compare. https://wisdem.readthedocs.io/en/master/examples/01_nrelcsm/tutorial.html
=#

using UnsteadyOpt
uo = UnsteadyOpt

Rtip = 63.0
D = 2 * Rtip
hub_height = 90.0
machine_rating = 5000.0
Vtip = 80.0
nblade = 3
n_bearing = 2

bos = uo.balance_of_station(machine_rating, D, hub_height)
#Wisdem example file has 2979e3, I've got 2019e3, sooo not good, but not terrible. 
bos_model = 2979e3 #Balance of station from the example file.

Q_rotor = uo.estimate_rotor_torque(machine_rating, D, Vtip)
m_blade = uo.estimate_blade_mass(D, 2, false) 
# m_blade = 26618.002029582505 #My current blade mass calculation. 

tcc = uo.calc_turbine_cost(machine_rating, D, Q_rotor, nblade, m_blade, n_bearing, hub_height)

icc = tcc + bos

# aep = 30660e6 # Example AEP value in kWh (from online calculator)
aep = 1.0187897726986462e10 #Current optimization original value in Wh #Todo. Should this be in kWh or MWh? -> Looks like it should just be in Wh. #Todo. Verify. -> Units, I think aep is suppoed to be in kWh, so my AEP value would be in the correct ballpark (maybe a little bit low). -> I think the example file has the cutout speed as 30 m/s... but I can't imagine that adds much. 
#The Wisdem example in the tutorial has 22290691.220265888. (2.229e7).... so potentially my aep is too high? 
aep_model = 22290691.220265888

o_m = uo.operating_maintenance_costs(aep)

lrc = uo.levelized_replacement_costs(machine_rating)

llc = uo.land_lease_costs(aep)

opex = o_m + llc + lrc #Example file has 144e3, which is significantly lower than the 90376.7392544637e3 from my model. 
opex_model = 144e3

tr = 0.4 #Tax rate
coe = uo.cost_of_energy(aep, bos, tcc, llc, lrc, o_m; tax_rate=tr)
coe_model = uo.cost_of_energy(aep_model, bos_model, tcc, opex_model; tax_rate=tr) #Note: This seems much more reasonable. If I just use the opex and the bos from the exampe, then figure out why my aep is so much higher.... then maybe I can just use these values. 

lcoe = uo.levelized_cost_of_energy(aep, bos, tcc, llc, lrc, o_m; discount_rate=0.07)

nothing