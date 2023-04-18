using BPM
using Test

# NOTE: The blunt trailing edge calculations have not yet been verified.  It would also
# be good to do a verification case for negative angles of attack.  Since this code
# is based on noise for a symmetric airfoil, the result should be symmetric.

@testset "Untripped Boundary Layer" begin

    # Verification against "Airfoil Self-Noise and Prediction", figure 43a

    # This verification covers untripped laminar and turbulent boundary layer noise

    # This test verifies our implementation by comparing the results of this code
    # with the results of the FORTRAN code written by Brooks, Pope, and Marcolini

    # frequencies
    f = [100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
        1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10000.0,
        12500.0, 16000.0, 20000.0, 25000.0, 31500.0, 40000.0]

    r_obs = 1.22       # observer distance
    θ_obs = 90*pi/180  # observer angles
    ϕ_obs = 90*pi/180  # observer angles
    L = 0.4572         # wetted span (m)
    c = 0.3048         # chord length (m)
    alpha = 1.516      # angle of attack (degrees)
    V = 71.3           # freestream velocity (m)
    c0 = 340.46        # air speed of sound (m/s)
    nu = 1.4529e-5     # air kinematic viscosity (m^2/s)
    tripped = false     # tripped or untripped boundary layer

    # temporary storage for pressure
    pressure = similar(f) .= 0

    # compute laminar boundary layer pressure contributions
    pressure .= 0
    BPM.add_laminar_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, L, c, alpha,
        V, c0, nu)
    spl_lam = @. 10*log10(pressure)

    # compute pressure side turbulent contributions
    pressure .= 0
    BPM.add_turbulent_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, L, c, alpha,
        V, c0, nu, tripped; pressure=true, suction=false, separation=false)
    spl_p = @. 10*log10(pressure)

    # compute suction side turbulent contributions
    pressure .= 0
    BPM.add_turbulent_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, L, c, alpha,
        V, c0, nu, tripped; pressure=false, suction=true, separation=false)
    spl_s = @. 10*log10(pressure)

    # compute separation turbulent contributions
    pressure .= 0
    BPM.add_turbulent_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, L, c, alpha,
        V, c0, nu, tripped; pressure=false, suction=false, separation=true)
    spl_a = @. 10*log10(pressure)

    # combine all contributions
    spl_tot = @. 10*log10(10.0^(spl_lam/10) + 10.0^(spl_p/10) + 10.0^(spl_s/10) + 10.0^(spl_a/10))

    # Computational Data from "Airfoil Self-Noise and Prediction", table D2

    data = [
        100.000 20.654 28.704 -100.000 -17.142 0.000 0.000 29.336
        125.000 24.461 31.965 -100.000 -13.285 0.000 0.000 32.676
        160.000 28.291 35.244 -75.254 -9.018 0.000 0.000 36.042
        200.000 31.437 37.937 -49.243 -5.161 0.000 0.000 38.815
        250.000 34.309 40.400 -27.506 -1.304 0.000 0.000 41.356
        315.000 37.023 42.736 -9.030 2.690 0.000 0.000 43.768
        400.000 39.577 44.949 6.266 6.820 0.000 0.000 46.057
        500.000 41.761 46.859 17.532 10.677 0.000 0.000 48.034
        630.000 43.845 48.706 26.603 14.671 0.000 0.000 49.954
        800.000 45.839 50.503 33.718 18.801 0.000 0.000 51.849
        1000.000 47.581 52.106 38.756 22.658 0.000 0.000 53.568
        1250.000 49.233 53.664 42.692 26.515 0.000 0.000 55.255
        1600.000 50.987 55.368 46.294 30.782 0.000 0.000 57.106
        2000.000 52.533 56.907 49.334 37.725 0.000 0.000 58.817
        2500.000 54.074 57.750 51.298 47.262 0.000 0.000 60.167
        3150.000 55.570 57.500 50.766 48.959 0.000 0.000 60.496
        4000.000 56.044 56.082 47.711 41.796 0.000 0.000 59.455
        5000.000 55.399 54.541 44.617 32.428 0.000 0.000 58.208
        6300.000 53.840 52.942 40.974 28.433 0.000 0.000 56.553
        8000.000 52.190 51.253 36.227 24.304 0.000 0.000 54.821
        10000.000 50.638 49.614 30.419 20.447 0.000 0.000 53.192
        12500.000 49.044 47.890 22.834 16.590 0.000 0.000 51.523
        16000.000 47.202 45.851 11.842 12.323 0.000 0.000 49.591
        20000.000 45.436 43.863 -0.924 8.466 0.000 0.000 47.731
        25000.000 43.549 41.710 -16.833 4.609 0.000 0.000 45.737
        31500.000 41.440 39.279 -37.092 0.614 0.000 0.000 43.503
        40000.000 39.065 36.522 -62.593 -3.515 0.000 0.000 40.987
    ]

    @test all(isapprox.(spl_lam, data[:,5], atol=0.25))
    @test all(isapprox.(spl_p, data[:,2], atol=0.25))
    @test all(isapprox.(spl_s, data[:,3], atol=0.25))
    @test all(isapprox.(max.(-100, spl_a), data[:,4], atol=0.25))
    @test all(isapprox.(spl_tot, data[:,8], atol=0.25))

    # # Uncomment to plot
    # using Plots
    # pyplot()

    # plot(minorticks=10, minorgrid=true,
    #     xaxis=:log,
    #     xlim = (0.2, 20.0),
    #     xticks = ([0.2, 1.0, 10.0, 20.0], [0.2, 1.0, 10.0, 20.0]),
    #     ylim = (40, 80)
    #     )

    # plot!(f/1000, spl_p, color=1, label="TBL-TE Pressure Side")
    # scatter!(data[:,1]/1000, data[:,2], color=1, label="")

    # plot!(f/1000, spl_s, color=2, label="TBL-TE Suction Side")
    # scatter!(data[:,1]/1000, data[:,3], color=2, label="")

    # plot!(f/1000, spl_a, color=3, label="Separation")
    # scatter!(data[:,1]/1000, data[:,4], color=3, label="")

    # plot!(f/1000, spl_tot, color=4, label="Total Prediction")
    # scatter!(data[:,1]/1000, data[:,8], color=4, label="")

    # plot!(f/1000, spl_lam, color=5, label="LBL-VS")
    # scatter!(data[:,1]/1000, data[:,5], color=5, label="")

end

@testset "Tripped Boundary Layer" begin

    # Verification against "Airfoil Self-Noise and Prediction", figure 91

    # This verification covers tripped boundary layer noise and tip noise

    # This test verifies our implementation by comparing the results of this code
    # with the results of the FORTRAN code written by Brooks, Pope, and Marcolini

    # frequency range
    f = [100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
    1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10000.0,
    12500.0, 16000.0, 20000.0, 25000.0, 31500.0, 40000.0]

    N = 10             # number of segments
    r_obs = 1.22       # observer distance
    θ_obs = 90*pi/180  # observer angles
    ϕ_obs = 90*pi/180  # observer angles
    L = 0.0305         # wetted span (m)
    c = 0.1524         # chord length (m)
    alpha = 5.4        # angle of attack (degrees)
    alpha_tip = 7.7    # tip angle of attack (degrees)
    V = 71.3           # freestream velocity (m)
    c0 = 340.46        # air speed of sound (m/s)
    nu = 1.4529e-5     # air kinematic viscosity (m^2/s)
    tripped = true     # tripped or untripped boundary layer
    round = true       # round wing tip
    aspect_ratio = 100 # aspect ratio

    # temporary storage for pressure
    pressure = similar(f) .= 0

    # compute tip vortex pressure contributions
    pressure .= 0
    BPM.add_tip_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, c, alpha_tip, V, c0, aspect_ratio, round)
    spl_tip = @. 10*log10(pressure)

    # compute pressure side turbulent contributions
    pressure .= 0
    BPM.add_turbulent_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, L, c, alpha,
        V, c0, nu, tripped; pressure=true, suction=false, separation=false)
    pressure *= N
    spl_p = @. 10*log10(pressure)

    # compute suction side turbulent contributions
    pressure .= 0
    BPM.add_turbulent_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, L, c, alpha,
        V, c0, nu, tripped; pressure=false, suction=true, separation=false)
    pressure *= N
    spl_s = @. 10*log10(pressure)

    # compute separation turbulent contributions
    pressure .= 0
    BPM.add_turbulent_pressure!(pressure, f, r_obs, θ_obs, ϕ_obs, L, c, alpha,
        V, c0, nu, tripped; pressure=false, suction=false, separation=true)
    pressure *= N
    spl_a = @. 10*log10(pressure)

    # combine contributions without tip
    spl_turb = @. 10*log10(10.0^(spl_p/10) + 10.0^(spl_s/10) + 10.0^(spl_a/10))

    # combine all contributions
    spl_tot = @. 10*log10(10.0^(spl_tip/10) + 10.0^(spl_p/10) + 10.0^(spl_s/10) + 10.0^(spl_a/10))

    # Computational Data from "Airfoil Self-Noise and Prediction", table D3

    data = [
        100.000 19.913 43.883 -19.803 0.000 0.000 -34.005 43.900
        125.000 23.788 46.159 -0.396 0.000 0.000 -24.312 46.184
        160.000 27.673 48.459 16.851 0.000 0.000 -14.255 48.498
        200.000 30.853 50.372 29.124 0.000 0.000 -5.769 50.452
        250.000 33.746 52.155 38.723 0.000 0.000 2.145 52.407
        315.000 36.470 53.894 46.334 0.000 0.000 9.738 54.662
        400.000 39.024 55.609 52.245 0.000 0.000 16.940 57.320
        500.000 41.202 57.165 56.460 0.000 0.000 23.074 59.897
        630.000 43.274 58.766 59.996 0.000 0.000 28.824 62.489
        800.000 45.252 60.360 63.297 0.000 0.000 34.121 65.130
        1000.000 46.980 60.940 65.719 0.000 0.000 38.475 67.016
        1250.000 48.620 60.473 65.697 0.000 0.000 42.257 66.917
        1600.000 50.364 58.874 62.909 0.000 0.000 45.774 64.582
        2000.000 51.911 57.328 59.818 0.000 0.000 48.349 62.363
        2500.000 53.456 55.775 56.383 0.000 0.000 50.351 60.580
        3150.000 54.709 54.122 51.975 0.000 0.000 51.821 59.364
        4000.000 54.799 52.336 45.974 0.000 0.000 52.694 58.443
        5000.000 53.761 50.565 38.550 0.000 0.000 52.917 57.439
        6300.000 52.162 48.597 28.510 0.000 0.000 52.544 56.204
        8000.000 50.507 46.387 15.081 0.000 0.000 51.512 54.736
        10000.000 48.936 44.132 -0.755 0.000 0.000 49.955 53.078
        12500.000 47.311 41.665 -20.241 0.000 0.000 47.826 51.110
        16000.000 45.415 38.655 -46.603 0.000 0.000 44.802 48.594
        20000.000 43.583 35.650 -75.275 0.000 0.000 41.466 46.075
        25000.000 41.611 32.347 -90.000 0.000 0.000 37.557 43.405
        31500.000 39.390 28.582 -90.000 0.000 0.000 32.904 40.555
        40000.000 36.873 24.291 -90.000 0.000 0.000 27.449 37.552
    ]

    @test all(isapprox.(spl_tip, data[:,7], atol=0.25))
    @test all(isapprox.(spl_p, data[:,2], atol=0.25))
    @test all(isapprox.(spl_s, data[:,3], atol=0.25))
    @test all(isapprox.(max.(-90, spl_a), data[:,4], atol=0.25))
    @test all(isapprox.(spl_tot, data[:,8], atol=0.25))

    # # Uncomment to plot
    # using Plots
    # pyplot()

    # plot(minorticks=10, minorgrid=true,
    #     xaxis=:log,
    #     xlim = (0.2, 20.0),
    #     xticks = ([0.2, 1.0, 10.0, 20.0], [0.2, 1.0, 10.0, 20.0]),
    #     ylim = (40, 90)
    #     )

    # plot!(f/1000, spl_tip, color=1, label="Tip")
    # scatter!(data[:,1]/1000, data[:,7], color=1, label="")

    # plot!(f/1000, spl_tot, color=2, label="Total with tip")
    # scatter!(data[:,1]/1000, data[:,8], color=2, label="")

    # plot!(f/1000, spl_turb, color=3, label="Total without tip")
    # scatter!(data[:,1]/1000, (@. 10*log10(10.0^(data[:,8]/10) - 10.0^(data[:,7]/10))), color=3, label="")

end

@testset "Rosiere Validation" begin

    # Rosiere Validation (243.84 m, 800 ft should be 47 dB)

    # This is a very rough test to see if the noise is reasonably close to the expected value.

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

    # Compute sound pressure levels
    oaspl, spl = BPM.sound_pressure_levels(ox, oy, oz, V, Ω, B, r, c, c1, h, alpha, psi, nu, c0;
        weighted=false, nbeta=1)

    println("Test Cases:")
    println("Rosiere Validation (47 dB): $oaspl")

    @test isapprox(oaspl, 47.0; atol=2.5)

end #Rosiere Test
