using FLOWMath, Plots
pyplot()

@testset "Smoothness" begin

    # NREL 5 MW Turbine
    data = [
        1.5    13.308  3.386 # Cylinder1.dat
        2.8667 13.308  3.542 # Cylinder1.dat
        5.6000 13.308  3.854 # Cylinder1.dat
        8.3333 13.308  4.167 # Cylinder2.dat
        11.7500 13.308 4.557 # DU40_A17.dat
        15.8500 11.480 4.652 # DU35_A17.dat
        19.9500 10.162 4.458 # DU35_A17.dat
        24.0500 9.011  4.249 # DU30_A17.dat
        28.1500 7.795  4.007 # DU25_A17.dat
        32.2500 6.544  3.748 # DU25_A17.dat
        36.3500 5.361  3.502 # DU21_A17.dat
        40.4500 4.188  3.256 # DU21_A17.dat
        44.5500 3.125  3.010 # NACA64_A17.dat
        48.6500 2.319  2.764 # NACA64_A17.dat
        52.7500 1.526  2.518 # NACA64_A17.dat
        56.1667 0.863  2.313 # NACA64_A17.dat
        58.9000 0.370  2.086 # NACA64_A17.dat
        61.6333 0.106  1.419 # NACA64_A17.dat
        63.0000 0.000  1.086 # NACA64_A17.dat
    ]

    # Scale rotor to 47 meters (to match the Vestas V47)
    r = data[:,1] * 47.0 / 126.0               # radial locations (m)
    c = data[:,3] * (47.0 / 126.0) * 2.190491  # chord length (m)
    c1 = 0.25*c                                # distance to leading edge from pitch axis (m)
    h = 0.01*c                                 # trailing edge thickness (m)
    alpha = data[:,2]                          # angle of attack (deg)
    psi = fill(14.0, length(r))                # solid angle (deg)

    # Define other inputs
    ox = 0.0    # x-location of observer relative to the turbine (m)
    oy = 243.84 # y-location of observer relative to the turbine (m)
    oz = -25.0    # z-location of observer relative to the turbine (m)
    V = 15.0    # wind velocity in the y-direction (m/s)
    Ω = 28.5 * 2*pi/60 # turbine rotation rate (rad/s)
    B = 3       # number of blades
    nu = 1.78e-5 # kinematic viscosity (m^2/s)
    c0 = 343.2 # speed of sound (m/s)

    # define function which modifies all the interesting variables
    Vfunc = (x) -> BPM.sound_pressure_levels(ox, oy, oz, x*V, Ω, B, r, c, c1, h, alpha, psi, nu, c0;
        turbulent=true, blunt=true, tip=false, weighted=false, smooth=false)[1]

    smooth_Vfunc = (x) -> BPM.sound_pressure_levels(ox, oy, oz, x*V, Ω, B, r, c, c1, h, alpha, psi, nu, c0;
        turbulent=true, blunt=true, tip=false, weighted=false, smooth=true)[1]

    x = range(2.0, 4.0, length=1000)
    plot(x*V, Vfunc.(x), label="",
        xlabel="Velocity (m/s)",
        ylabel="OASPL")
    plot!(x*V, smooth_Vfunc.(x), label="")

    afunc = (x) -> BPM.sound_pressure_levels(ox, oy, oz, V, Ω, B, r, c, c1, h, x.*alpha, psi, nu, c0;
        turbulent=true, blunt=true, tip=false, weighted=false, smooth=false)[1]

    smooth_afunc = (x) -> BPM.sound_pressure_levels(ox, oy, oz, V, Ω, B, r, c, c1, h, x.*alpha, psi, nu, c0;
        turbulent=true, blunt=true, tip=false, weighted=false, smooth=true)[1]

    x = range(0.0, 2.0, length=1000)
    plot(x * alpha[1], afunc.(x), label="",
        xlabel="Root Angle of Attack (deg)",
        ylabel="OASPL")
    plot!(x * alpha[1], smooth_afunc.(x), label="")

end #Rosiere Test
