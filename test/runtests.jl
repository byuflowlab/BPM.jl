using BPM
using Base.Test

@testset "Rosiere Validation" begin
##################################################################################
##################################################################################
##################################################################################

# Rosiere Valdiation (243.84 m, 800 ft should be 47 dB)
# http://www.mge.com/environment/green-power/wind/rosiere.htm
x_test = [0.] # x-locations of turbines (m)
y_test = [0.] # y-locations of turbines (m)
obs_test = [0., 243.84, 0.] # x-, y-, and z-location of the observer (m)
winddir_test = 0. # wind direction (deg)
rpm_test = [28.5] # rotation rate of the tubrines (rpm)
windvel_test = [15.] # wind velocity (m/s)
B_test = 3 # number of blades
h_test = 25. # height of the turbine hub (m)
noise_corr = 0.8697933840957954 # correction factor for noise

rad =  [1.069324603174603, 2.088888888888889, 3.1084531746031745, 4.382936507936508, 5.912301587301587, 7.441666666666666, 8.971031746031747, 10.500396825396825, 12.029761904761905, 13.559126984126985, 15.088492063492065, 16.617857142857144, 18.147222222222222, 19.6765873015873, 20.951070634920637, 21.97063492063492, 22.990199206349207, 23.5] # radial positions (m)
c =  [2.8941253867777776, 3.1490568155396828, 3.404805332214286, 3.7234696181666673, 3.8010929698730163, 3.6425779148095243, 3.4718065410555554, 3.2740712661825397, 3.062445496793651, 2.861441870269841, 2.660438243746032, 2.459434617222222, 2.258430990698413, 2.057427364174603, 1.889924342071429, 1.7044453858888888, 1.1594477481190477] # chord lengths (m)
alpha =  [13.308, 13.308, 13.308, 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.37, 0.106] # angles of attack (deg)
c1 = c*0.25 # pitch axis (m)
AR = 17. # blade aspect ratio

nu = 1.78e-5 # kinematic viscosity (m^2/s)
c0 = 343.2 # speed of sound (m/s)
psi = 14.0 # solid angle (deg)
# f    (nturb,  nseg,     nobs,       x,            y,          obs,    winddir,windvel,rpm, B, Hub, rad,  c,  c1,alpha,nu,    c0,      psi, AR, noise_corr)
db_test_ros = BPM.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, B_test, h_test, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)

println("Test Cases:")
println("Rosiere Validation (47 dB): $db_test_ros")

@test isapprox(db_test_ros, 47.0; atol=1e-6)

end #Rosiere Test


@testset "Eric's Example" begin
##################################################################################
##################################################################################
##################################################################################

# SPL Test Paremeters (changed to whatever desired)
x_test = [0.,5.] # x-locations of turbines (m)
y_test = [0.,5.] # y-locations of turbines (m)
obs_test = [0., 200., 0.] # x-, y-, and z-location of the observer (m)
winddir_test = 0. # wind direction (deg)
rpm_test = [28.5,28.5] # rotation rate of the tubrines (rpm)
windvel_test = [10.0,10.0] # wind velocity (m/s)
B_test = 3 # number of blades
h_test = 25. # height of the turbine hub (m)
noise_corr = 0.8697933840957954 # correction factor for noise

rad =  [1.069324603174603, 2.088888888888889, 3.1084531746031745, 4.382936507936508, 5.912301587301587, 7.441666666666666, 8.971031746031747, 10.500396825396825, 12.029761904761905, 13.559126984126985, 15.088492063492065, 16.617857142857144, 18.147222222222222, 19.6765873015873, 20.951070634920637, 21.97063492063492, 22.990199206349207, 23.5] # radial positions (m)
c =  [2.8941253867777776, 3.1490568155396828, 3.404805332214286, 3.7234696181666673, 3.8010929698730163, 3.6425779148095243, 3.4718065410555554, 3.2740712661825397, 3.062445496793651, 2.861441870269841, 2.660438243746032, 2.459434617222222, 2.258430990698413, 2.057427364174603, 1.889924342071429, 1.7044453858888888, 1.1594477481190477] # chord lengths (m)
alpha =  [13.308, 13.308, 13.308, 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.37, 0.106] # angles of attack (deg)
c1 = c*0.25 # pitch axis (m)
AR = 17. # blade aspect ratio

nu = 1.78e-5 # kinematic viscosity (m^2/s)
c0 = 343.2 # speed of sound (m/s)
psi = 14.0 # solid angle (deg)

db_test = BPM.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, B_test, h_test, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)

println("Test SPL (50.8446366094): $db_test")
@test isapprox(db_test, 50.8446366094; atol=1e-6)

end #Eric's other example

@testset "Mariah Windspire 5 MW" begin
##################################################################################
##################################################################################
##################################################################################

div=5

nu = 1.78e-5 # kinematic viscosity (m^2/s)
c0 = 343.2 # speed of sound (m/s)
psi = 14.0 # solid angle (deg)

turbx = [0.] # turbine x-positions (m)
turby = [0.] # turbine y-positions (m)
Vinf = 8. # free stream wind speed (m/s)
dia = 1.2 # turbine diameter (m)
rad = dia/2. # turbine radius (m)
tsrd = 2.625 # tip-speed ratio

rot = [1]*tsrd*Vinf/rad # turbine rotation rates (rad/s)

winddir = 180. # wind direction (deg)

# Aerodynamic properties
twist = 0.0
delta = 0.0

B = 3 # number of blades
chord = 0.128 # chord length (m)
c = ones(div)*chord # chord lengths over height positions (m)
c1 = c*0.5 # pitch axis location (m)
alpha = ones(div)*0.0 # angles of attack (deg)
Hub = 2. # hub height (m)
H = 6.1 # blade height (m)
high = linspace(0,H,div+1) # height positions of the blade (m)

AR = 5. # aspect ratio

# Wake velocities from surronding turbines
ntheta = 4
wakex = zeros(ntheta)
wakey = zeros(ntheta)

noise_corr = 1.0 # correction factor for noise

#Test Point
X=5
Y=5
Z=0
SPL=BPM.turbinepos_VAWT(ntheta, turbx, turby, [X,Y,Z], winddir, B, Hub, high, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr, rot, Vinf, wakex, wakey)

println("Test SPL (70.3): $SPL")
@test isapprox(db_test, 70.263597969; atol=1e-6)

end
