using Test, UnsteadyOpt, CCBlade

uo = UnsteadyOpt


@testset "aerodynamics.jl" begin

    airfoilpath = "/Users/adamcardoza/repos/Cardoza2024_Unsteady_Aerostructural_Blade_Opt/src/optimization/static_aerostructural/clt/airfoils_interpolated"

    # airfoil = uo.allocate_interp_airfoil(1, airfoilpath)
    # @code_warntype uo.allocate_interp_airfoil(1, airfoilpath) #Note: Because of the readdlm, this is inherently type unstable.
    # @show @inferred(uo.allocate_interp_airfoil(1, airfoilpath)) #Todo: Erroring. 
    # @test isa(@inferred(uo.allocate_interp_airfoil(1, airfoilpath)), CCBlade.Cylinder)

    rvec = [1.80306613137, 1.90155987557, 2.00005361977, 2.10470322299, 2.20319696719, 2.30169071139, 2.86802974054, 3.00345863882, 3.10195238302, 5.60738700113, 7.00476699698, 8.3405884027, 10.5074507751, 11.7632460137, 13.5115099732, 15.863048116, 18.5162233505, 19.9690060774, 22.0189071286, 24.0749640388, 26.12486509, 28.1747661411, 32.2807241025, 33.5303634821, 36.3866820639, 38.5350768593, 40.4864841663, 42.5425410764, 43.5397902365, 44.5924421276, 46.5438494346, 48.698400089, 52.7982021914, 56.2208598023, 58.9540612039, 61.6934184645]
    airfoils = uo.get_interp_polars(rvec, airfoilpath)

    Vmean = 6.0
    Vin = 3.0
    Vout = 25.0
    Vcurve = collect(Vin:1.0:Vout)
    nwind = length(Vcurve)

    Rhub = 1.5
    Rtip = 63.0
    B = 3
    rotor = Rotor(Rhub, Rtip, B, turbine=true)

    tsr0 = 7.0175
    pitch = 0.0
    pitches = ones(nwind) .* pitch

    yaw = 0.0
    tilt = 0.0
    azimuth = 0.0
    hubHt = 90.0
    shearExp = 0.2
    rho = 1.225

    cvec = [3.30974143717, 3.32552109037, 3.34126379956, 3.35794850262, 3.37361093258, 3.38923246006, 3.47819767927, 3.49923810816, 3.51447938721, 3.88091952034, 4.06260542204, 4.21678354504, 4.41794897837, 4.50221292399, 4.57390637252, 4.57496258253, 4.49355793249, 4.44296484557, 4.36449116824, 4.27745048448, 4.18236399597, 4.07898355526, 3.84696806115, 3.7697498699, 3.58167705636, 3.42960605033, 3.28382135203, 3.12880081291, 3.05388779827, 2.97485038083, 2.82801258654, 2.66459194226, 2.34610078292, 2.06819798348, 1.83533535858, 1.58981373477]
    twistvec = [0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.23173655953085248, 0.22036544910113506, 0.20545837279536666, 0.18915019829733667, 0.18044007588271704, 0.1684029212782924, 0.1566149422632624, 0.14513335086769033, 0.13390909268418555, 0.11231246579079515, 0.10601707527729945, 0.09211716927640182, 0.08211110735539948, 0.0733564662404205, 0.06447627739537955, 0.060296223298632066, 0.05597401823528258, 0.04821186406774127, 0.040118225037764965, 0.026132878944916376, 0.015797360903874846, 0.008362959903121539, 0.0015980174008993844]
    sections = CCBlade.Section.(rvec, cvec, twistvec, airfoils)

    TF = Float64

    pc, aep = uo.get_powercurve(rotor, sections, Vcurve, tsr0, pitches, yaw, tilt, azimuth, hubHt, shearExp, rho, Vmean, TF)

    # @code_warntype uo.get_powercurve(rotor, sections, Vcurve, tsr0, pitches, yaw, tilt, azimuth, hubHt, shearExp, rho, Vmean, TF) #todo: potentially type unstable because of the rvec and the solve... 
    @test isa(@inferred(uo.get_powercurve(rotor, sections, Vcurve, tsr0, pitches, yaw, tilt, azimuth, hubHt, shearExp, rho, Vmean, TF)), Tuple{Vector{Float64}, Float64})

end