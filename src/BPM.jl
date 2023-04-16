module BPM

using FLOWMath

export sound_pressure_levels

# Default one-third octave band frequencies (Hz)
const default_f = [  100.0,   125.0,   160.0,   200.0,   250.0,   315.0,   400.0,   500.0,
                     630.0,   800.0,  1000.0,  1250.0,  1600.0,  2000.0,  2500.0,  3150.0,
                    4000.0,  5000.0,  6300.0,  8000.0, 10000.0, 12500.0, 16000.0, 20000.0,
                   25000.0, 31500.0, 40000.0]

# A-weighting curve (dBA) for sound perception correction
const default_AdB = [-19.145, -16.190, -13.244, -10.847, -8.675, -6.644, -4.774, -3.248,
                      -1.908,  -0.795,   0.0,     0.576,  0.993,  1.202,  1.271,  1.202,
                       0.964,   0.556,  -0.114,  -1.144, -2.488, -4.250, -6.701, -9.341,
                     -12.322, -15.694, -19.402]

"""
    sound_pressure_levels(ox, oy, oz, V, Ω, B, r, c, c1, h, alpha, psi, nu, c0; kwargs...)

Calculate the (A-weighted) sound pressure levels (SPL) of a turbine.  Return the overall
sound pressure level (OASPL) and a vector containing sound pressure levels (SPL) for each
frequency.

# Arguments
 - `ox`: Observer x-location (relative to the turbine)
 - `oy`: Observer y-location (relative to the turbine)
 - `oz`: Observer z-location (relative to the turbine)
 - `V`: Wind velocity (in the y-direction) (m/s)
 - `Ω`: Rotation rate (rad/s)
 - `B`: Number of blades
 - `r`: Vector of blade radial locations (m)
 - `c`: Chord lengths at each radial location (m)
 - `c1`: Distance from the pitch axis to the leading edge at each radial location (m)
 - `h`: Trailing edge thickness at each radial location (m)
 - `alpha`: Angle of attack at each radial location (degrees)
 - `psi`: Wedge angle at each radial location (degrees)
 - `nu`: air kinematic viscosity (m^2/s)
 - `c0`: air speed of sound (m/s)

# Keyword Arguments
 - `laminar = false`: Flag(s) to compute laminar boundary layer noise
 - `turbulent = true`: Flag(s) to compute turbulent boundary layer noise
 - `blunt = true`: Flag(s) to compute blunt trailing edge noise
 - `tip = true`: Flag(s) to compute tip noise
 - `trip = true`: Flag(s) to trip boundary layer
 - `round = true`: Flag which indicates whether the tip is round or flat
 - `weighted = false`: Indicates whether to apply an A-weighting to the sound pressure levels
 - `nbeta = ceil(Int, 8/B)`: Number of rotation angles to consider.
 - `f = BPM.default_f`: Frequencies to use in the analysis
 - `Adb = BPM.default_AdB`: A-weighting for each frequency
"""
function sound_pressure_levels(ox, oy, oz, V, Ω, B, r, c, c1, h, alpha, psi, nu, c0;
    laminar=false, turbulent=true, blunt=true, tip=true, trip=true, round=true, weighted=true,
    nbeta = ceil(Int, 8/B), f = default_f, AdB = default_AdB)

    # get floating point type
    TF = promote_type(typeof(ox), typeof(oy), typeof(oz), typeof(V), typeof(Ω), eltype(r),
        eltype(c), eltype(c1), eltype(alpha), eltype(psi), typeof(nu), typeof(c0))

    # check input lengths
    @assert length(r) == length(c) == length(c1) == length(alpha) == length(psi)
    @assert !weighted || length(f) == length(default_AdB)

    # number of radial stations
    nr = length(r)
    nf = length(f)

    # define blade angles
    beta = range(0.0, 2*pi/B, nbeta+1)[1:end-1]

    # calculate aspect ratio
    area = FLOWMath.trapz(r, c)
    cbar = area/(r[end]-r[1])
    aspect_ratio = r[end]/cbar

    # calculate tip velocity
    Vt = sqrt(V^2 + (Ω*r[end])^2)

    # initialize sound pressure levels for each frequency
    spl = zeros(TF, nf)

    # initialize pressure for each frequency
    p = zeros(TF, nf)

    # loop through each rotation increment
    for ibeta = 1 : nbeta

        # reset pressure accumulation array
        p .= 0

        # loop through each blade
        for iB = 1 : B

            # current blade angle
            cbeta = beta[ibeta] + (iB - 1) * 2*pi/B

            if tip
                # find observer location relative to the tip
                r_obs, θ_obs, ϕ_obs = observer_location(ox, oy, oz, c[end], c1[end], r[end], cbeta)
                # add tip vortex pressure contributions
                add_tip_pressure!(p, f, r_obs, θ_obs, ϕ_obs, c[end], alpha[end], Vt, c0,
                    aspect_ratio, round)
            end

            # loop through each radial section
            for k = 1:nr-1

                # calculate blade section properties
                Lm = r[k+1] - r[k]             # segment length (m)
                rm = (r[k+1] + r[k])/2         # segment radius (m)
                cm = (c[k+1] + c[k])/2         # segment chord (m)
                c1m = (c1[k+1] + c1[k])/2      # distance to segment leading edge (m)
                hm = (h[k+1] + h[k])/2         # segment trailing edge thickness (m)
                am = (alpha[k+1] + alpha[k])/2 # segment angle of attack (deg)
                pm = (psi[k+1] + psi[k])/2     # segment wedge angle (deg)
                Vm = sqrt(V^2 + (Ω*rm)^2)      # segment local velocity (m/s)

                # find observer location relative to this blade element
                r_obs, θ_obs, ϕ_obs = observer_location(ox, oy, oz, cm, c1m, rm, cbeta)

                # check if boundary layer is tripped on this airfoil
                tripped = typeof(trip) <: AbstractVector ? trip[k] : trip

                if typeof(laminar) <: AbstractVector ? laminar[k] : laminar
                    # add laminar boundary layer pressure contributions
                    add_laminar_pressure!(p, f, r_obs, θ_obs, ϕ_obs, Lm, cm, am, Vm, c0, nu, tripped)
                end

                if typeof(turbulent) <: AbstractVector ? turbulent[k] : turbulent
                    # add turbulent boundary layer pressure contributions
                    add_turbulent_pressure!(p, f, r_obs, θ_obs, ϕ_obs, Lm, cm, am, Vm, c0, nu, tripped)
                end

                if typeof(blunt) <: AbstractVector ? blunt[k] : blunt
                    # add blunt trailing edge pressure contributions
                    add_bluntness_pressure!(p, f, r_obs, θ_obs, ϕ_obs, Lm, cm, hm, am, pm, Vm, c0, nu, tripped)
                end

            end
        end

        # accumulate squared pressures at each frequency
        @. spl += p^2

    end

    # compute root mean squared pressure
    @. spl = sqrt(spl/nbeta)

    # convert to decibels
    @. spl = 10.0*log10(spl)

    # apply A-weighting
    if weighted
        spl .+= AdB
    end

    # compute overall sound pressure level
    psum2 = zero(TF)
    for i = 1:length(spl)
        psum2 += (10.0 ^ (2*spl[i]/10.0))
    end
    oaspl = 10*log10(sqrt(psum2))

    # return noise in decibels at each frequency
    return oaspl, spl
end

# computed directivity angles and distance for a HAWT
# based on work by Luis Vargas (Wind Turbine Noise Prediction)
function observer_location(xo, yo, zo, c, c1, d, beta)

    # Calculate the trailing edge position relative to the hub
    xs = sin(beta)*d - cos(beta)*(c - c1)
    zs = cos(beta)*d + sin(beta)*(c - c1)

    # Calculate the observer position relative to the trailing edge
    xe_d = xo - xs
    ze_d = zo - zs

    # Rotate the observer position by the blade angle (beta)
    theta = pi - beta
    xe =  cos(theta)*xe_d + sin(theta)*ze_d
    ze = -sin(theta)*xe_d + cos(theta)*ze_d

    # Calculate the distance to the observer and the directivity angles
    r = sqrt(yo^2+xe^2+ze^2)
    θ = atan(sqrt(yo^2+ze^2), xe)
    ϕ = atan(yo, ze)

    # Quadratic smoothing when ϕ is close to 0 or 180 degrees
    if (abs(ϕ) < 5.0*pi/180.0)
        if ϕ >= 0.0
            sign = 1
        else
            sign = -1
        end
        phi_er = abs(ϕ)*180.0/pi
        phi_er = 0.1*phi_er^2+2.5
        ϕ = sign*phi_er*pi/180.0
    elseif (abs(ϕ) > 175.0*pi/180.0)
        if ϕ >= 0.0
            sign = 1
        else
            sign = -1
        end
        phi_er = abs(ϕ)*180.0/pi
        phi_er = -0.1*(phi_er-180.0)^2+177.5
        ϕ = sign*phi_er*pi/180.0
    end

    return r, θ, ϕ
end

# Laminar Boundary Layer Vortex Shedding Noise (pg. 66)
function add_laminar_pressure!(p, f, r, θ, ϕ, L, c, alpha, V, c0, nu, trip)

    # Mach Number
    M = V/c0

    # Reynolds Number (based on chord)
    Re = V*c/nu

    # compute boundary layer thickness
    dp = boundary_thickness(Re, c, alpha)

    # compute directivity function
    Dh = Dhfunc(M, θ, ϕ)

    # compute reference Strouhal number (eq. 54)
    if Re <= 1.3e5
        St1p = 0.18
    elseif 1.3e5 < Re <= 4.0e5
        St1p = 0.001756*Re^0.3931
    else # 4.0e5 < Re
        St1p = 0.28
    end

    # compute peak Strouhal (eq. 56)
    Stp_peak = St1p*10.0^(-0.04*alpha)

    # compute reference Reynolds number
    if alpha <= 3.0
        Re0 = 10.0^(0.215*alpha + 4.978)
    else # 3.0 < alpha
        Re0 = 10.0^(0.120*alpha + 5.263)
    end

    # compute peak scaled spectrum level
    d = Re/Re0

    if d <= 0.3237
        G2 = 77.852*log10(d) + 15.328
    elseif 0.3237 < d <= 0.5689
        G2 = 65.188*log10(d) + 9.125
    elseif 0.5689 < d <= 1.7579
        G2 = -114.052*log10(d)^2
    elseif 1.7579 < d <= 3.0889
        G2 = -65.188*log10(d) + 9.125
    else # 3.0889 < d
        G2 = -77.852*log10(d) + 15.328
    end

    G3 = 171.04 - 3.03*alpha

    scale = 10.0*log10(dp*M^5*L*Dh/r^2)

    # compute pressure levels for each frequency
    for i = 1:length(f)

        Stp = f[i]*dp/V

        e = Stp/Stp_peak

        if e <= 0.5974
            G1 = 39.8*log10(e) - 11.12
        elseif 0.5974 < e <= 0.8545
            G1 = 98.409*log10(e) + 2.0
        elseif 0.8545 < e <= 1.17
            G1 = -5.076 + sqrt(2.484 - 506.25*log10(e)^2)
        elseif 1.17 < e <= 1.674
            G1 = -98.409*log10(e) + 2.0
        else # 1.674 < e
            G1 = -39.8*log10(e) - 11.12
        end

        spl_lam = G1 + G2 + G3 + scale

        p[i] += 10.0^(spl_lam/10.0)
    end

    return p
end

# Turbulent Boundary Layer Trailing Edge Noise and Separated Flow Noise (pg. 59)
function add_turbulent_pressure!(p, f, r, θ, ϕ, L, c, alpha, V, c0, nu, trip;
    pressure=true, suction=true, separation=true)

    # Mach number
    M = V/c0

    # Reynolds number (based on chord)
    Re = V*c/nu

    # compute boundary layer thickness
    dp, ds = displacement_thickness(Re, c, alpha, trip)

    # directivity functions
    Dl = Dlfunc(M, θ, ϕ)
    Dh = Dhfunc(M, θ, ϕ)

    # compute Reynolds numbers based on displacement thicknesses
    Rp = V*dp/nu
    Rs = V*ds/nu

    # determine peak Strouhal numbers for `A` and `B` curve calculations
    St1 = 0.02*M^(-0.6)

    if alpha <= 1.333
        St2 = St1
    elseif 1.333 < alpha <= 12.5
        St2 = St1*10.0^(0.0054*(alpha - 1.333)^2)
    else # 12.5 < alpha
        St2 = St1*4.72
    end

    St1bar = (St1 + St2)/2.0

    # construct A-curve
    a0 = a0_func(Re)
    Amin0 = Amin_func(a0)
    Amax0 = Amax_func(a0)
    A_ratio = (20 + Amin0)/(Amin0 - Amax0)

    # construct alternative A-curve
    ap0 = a0_func(3*Re)
    Apmin0 = Amin_func(ap0)
    Apmax0 = Amax_func(ap0)
    Ap_ratio = (20 + Apmin0)/(Apmin0 - Apmax0)

    # construct B-curve
    b0 = b0_func(Re)
    Bmin0 = Bmin_func(b0)
    Bmax0 = Bmax_func(b0)
    B_ratio = (20.0 + Bmin0)/(Bmin0 - Bmax0)

    # compute K1
    if Re < 2.47e5
        K1 = -4.31*log10(Re) + 156.3
    elseif 2.47e5 <= Re < 8.0e5
        K1 = -9.0*log10(Re) + 181.6
    else # 8.0e5 < Re
        K1 = 128.5
    end

    # compute ΔK1 (eq. 48)
    if Rp <= 5000.0
        ΔK1 = -alpha*(5.29 - 1.43*log10(Rp))
    else
        ΔK1 = 0.0
    end

    # compute K2 (eq. 49)
    gamma =  27.094*M +  3.31
    beta =   72.650*M + 10.74
    gamma0 = 23.430*M +  4.651
    beta0 = -34.190*M - 13.820

    if alpha <= gamma0 - gamma
        K2 = K1 - 1000.0
    elseif gamma0 - gamma < alpha <= gamma0 + gamma
        K2 = K1 + sqrt(beta^2 - (beta/gamma)^2*(alpha-gamma0)^2) + beta0
    else
        K2 = K1 - 12.0
    end

    # check if the airfoil is stalling/stalled
    if alpha > 12.5 || alpha >= gamma0

        # compute pressure for each frequency
        for i = 1:length(f)

            # pressure side (not included in theory document)
            spl_p = 10.0*log10(dp*M^5*L*Dl/r^2)

            # suction side (not included in theory document)
            spl_s = 10.0*log10(ds*M^5*L*Dl/r^2)

            # angle of attack
            Sts = f[i]*ds/V
            b = abs(log10(Sts/St2))
            Amin = Amin_func(b)
            Amax = Amax_func(b)
            Aa = Amin + Ap_ratio*(Amax - Amin)
            spl_a = Aa + K2 + 10.0*log10(ds*M^5*L*Dl/r^2)

            # add to pressure for each frequency
            pressure && (p[i] += 10^(spl_p/10))
            suction && (p[i] += 10^(spl_s/10))
            separation && (p[i] += 10^(spl_a/10))
        end

    else

        # compute pressure for each frequency
        for i = 1:length(f)

            # pressure side
            Stp = f[i]*dp/V
            ap = abs(log10(Stp/St1))
            Apmin = Amin_func(ap)
            Apmax = Amax_func(ap)
            Ap = Apmin + A_ratio*(Apmax - Apmin)
            spl_p = 10.0*log10(dp*M^5*L*Dh/r^2) + Ap + K1 - 3 + ΔK1

            # suction side
            Sts = f[i]*ds/V
            as = abs(log10(Sts/St1bar)) # abs(log10(Sts/St1))
            Asmin = Amin_func(as)
            Asmax = Amax_func(as)
            As = Asmin + A_ratio*(Asmax - Asmin)
            spl_s = 10.0*log10(ds*M^5*L*Dh/r^2) + As + K1 - 3

            # angle of attack
            b = abs(log10(Sts/St2))
            Bmin = Bmin_func(b)
            Bmax = Bmax_func(b)
            Ba = Bmin + B_ratio*(Bmax - Bmin)
            spl_a = 10.0*log10(ds*M^5*L*Dh/r^2) + Ba + K2

            # add to pressure for each frequency
            pressure && (p[i] += 10^(spl_p/10))
            suction && (p[i] += 10^(spl_s/10))
            separation && (p[i] += 10^(spl_a/10))

        end

    end

    return p
end

"""
    Amin_func(a)

Defines the curve fit corresponding to the A-curve for the minimum allowed Reynolds number
"""
function Amin_func(a)

    a = abs(a)

    if a < 0.204
        Amin = sqrt(67.552 - 886.788*a^2) - 8.219
    elseif 0.204 <= a <= 0.244
        Amin = -32.665*a + 3.981
    else # 0.244 > a
        Amin = -142.795*a^3 + 103.656*a^2-57.757*a+6.006
    end

    return Amin
end

"""
    Amax_func(a)

Defines the curve fit corresponding to the A-curve for the maximum allowed Reynolds number
"""
function Amax_func(a)

    a = abs(a)

    if a < 0.13
        Amax = sqrt(67.552 - 886.788*a^2) - 8.219
    elseif 0.13 <= a <= 0.321
        Amax = -15.901*a + 1.098
    else # 0.321 < a
        Amax = -4.669*a^3 + 3.491*a^2 - 16.699*a + 1.149
    end

    return Amax
end

"""
    Bmin_func(a)

Defines the curve fit corresponding to the B-curve for the minimum allowed Reynolds number
"""
function Bmin_func(b)

    b = abs(b)

    if b < 0.13
        Bmin = sqrt(16.888 - 886.788*b^2) - 4.109
    elseif 0.13 <= b <= 0.145
        Bmin = -83.607*b + 8.138
    else # 0.145 < b
        Bmin = -817.810*b^3 + 355.210*b^2 - 135.024*b + 10.619
    end

    return Bmin
end

"""
    Bmax_func(a)

Defines the curve fit corresponding to the B-curve for the maximum allowed Reynolds number
"""
function Bmax_func(b)

    b = abs(b)

    if b < 0.10
        Bmax = sqrt(16.888 - 886.788*b^2) - 4.109
    elseif 0.10 <= b <= 0.187
        # Bmax = -31.330*b + 1.854 # <-- version in theory
        Bmax = -31.313*b + 1.854 # <-- version in code
    else # 0.187 < b
        Bmax = -80.541*b^3 + 44.174*b^2 - 39.381*b + 2.344
    end

    return Bmax
end

"""
    a0_func(Re)

Determines where the A-curve takes on a value of -20 dB
"""
function a0_func(Re)

    if Re < 9.52e4
        a0 = 0.57
    elseif 9.52e4 <= Re <= 8.57e5
        a0 = -9.57e-13*(Re - 8.57e5)^2 + 1.13
    else # 8.57e5 < Re
        a0 = 1.13
    end

    return a0
end

"""
    b0_func(Re)

Determines where the B-curve takes on a value of -20 dB
"""
function b0_func(Re)

    if Re < 9.52e4
        b0 = 0.30
    elseif 9.52e4 <= Re < 8.57e5
        b0 = -4.48e-13*(Re - 8.57e5)^2 + 0.56
    else # 8.57e5 < Re
        b0 = 0.56
    end

    return b0
end

"""
    Dhfunc(M, θ, ϕ)

Computes the high frequency directivity function for the input observer location.
"""
function Dhfunc(M, θ, ϕ)

    # NOTE: This function should not be used for high-angle separation since it becomes
    # inaccurate as θ approaches 180 deg

    # assumed convection Mach number (pg. 107)
    Mc = 0.8*M

    # directivity function (eq. B1)
    Dh = 2*sin(θ/2.0)^2*sin(ϕ)^2/((1 + M*cos(θ)) * (1 + (M - Mc)*cos(θ))^2)

    return Dh
end

"""
    Dlfunc(M, θ, ϕ)

Computes the low-frequency directivity function for the input observer location.
"""
function Dlfunc(M, θ, ϕ)

    # directivity function (eq. B2)
    Dl = sin(θ)^2*sin(ϕ)^2/(1 + M*cos(θ))^4

    return Dl
end

# Trailing Edge Bluntness Vortex Shedding Noise (pg. 78)
function add_bluntness_pressure!(p, f, r, θ, ϕ, L, c, h, alpha, psi, V, c0, nu, trip)

    # Mach number
    M = V/c0

    # Reynolds number
    Re = (V*c)/nu

    # compute displacement thickness
    dp, ds = displacement_thickness(Re, c, alpha, trip)

    # compute average displacement thickness (eq. 73)
    dav = (dp + ds)/2.0

    # ratio of trailing edge thickness to average boundary-layer displacement thickness
    hdav = h/dav

    # directivity function (eq. B1)
    Dh = Dhfunc(M, θ, ϕ)

    # peak Strouhal number (eq 72)
    if 0.2 <= hdav
        St_peak = (0.212 - 0.0045*psi)/(1.0 + 0.235/hdav - 0.0132/hdav^2)
    else # hdav < 0.2
        St_peak = 0.1*hdav + 0.095 - 0.00243*psi
    end

    # compute scaled spectrum level
    if hdav <= 5
        G4 = 17.5*log10(hdav) + 157.5 - 1.114*psi
    else # 5 < hdav
        G4 = 169.7 - 1.114*psi
    end

    # compute pressure for each frequency
    for i = 1:length(f)

        St = (f[i]*h)/V

        eta = log10(St/St_peak)

        # calculate G5 at psi = 0 deg

        hdav_prime = 6.724*hdav^2 - 4.019*hdav + 1.107

        G5_0 = G5func(hdav_prime, eta)

        # calculate G5 at psi = 14 deg

        G5_14 = G5func(hdav, eta)

        # linearly interpolate between 0 degrees and 14 degrees (equation 75)
        G5 = G5_0 + 0.0714*psi*(G5_14 - G5_0)

        # this check is not found in the theory, but exists in the code
        tmp = G5func(0.25, eta)
        if G5 > 0
            G5 = 0
        elseif G5 > tmp
            G5 = tmp
        end

        scale = 10.0*log10((M^5.5*h*Dh*L)/r^2)

        spl_blunt = G4 + G5 + scale

        p[i] += 10^(spl_blunt/10)

    end

    return p
end

function G5func(hdav, eta)

    # equation 78
    if hdav < 0.25
        mu = 0.1221
    elseif 0.25 <= hdav < 0.62
        mu = -0.2175*hdav + 0.1755
    elseif 0.62 <= hdav < 1.15
        mu = -0.0308*hdav + 0.0596
    else # 1.15 <= hdav
        mu = 0.0242
    end

    # equation 79
    if hdav <= 0.02
        m = 0.0
    elseif 0.02 < hdav <= 0.5
        m = 68.724*hdav - 1.35
    elseif 0.5 < hdav <= 0.62
        m = 308.475*hdav - 121.23
    elseif 0.62 < hdav <= 1.15
        m = 224.811*hdav - 69.35
    elseif 1.15 < hdav <= 1.2
        m = 1583.28*hdav - 1631.59
    else # 1.2 < hdav
        m = 268.344
    end

    # equation 80
    eta_0 = -sqrt((m^2*mu^4)/(6.25 + m^2*mu^2))

    # equation 81
    k = 2.5*sqrt(1.0 - (eta_0/mu)^2) - 2.5 - m*eta_0

    # equation 76
    if eta < eta_0
        G5 = m*eta+k
    elseif eta_0 <= eta < 0.0
        G5 = 2.5*sqrt(1.0 - (eta/mu)^2) - 2.5
    elseif 0.0 <= eta < 0.03616
        G5 = sqrt(1.5625 - 1194.99*eta^2) - 1.25
    else # 0.03616 <= eta
        G5 = -155.543*eta + 4.375
    end

    return G5
end

# tip vortex noise correction data based on "Airfoil Tip Vortex Formation Noise"
const aspect_data = [2.0,2.67,4.0,6.0,12.0,24.0]
const aratio_data = [0.54,0.62,0.71,0.79,0.89,0.95]
const aspect_ratio_correction = FLOWMath.Akima(aspect_data, aratio_data)

# Tip Vortex Formation Noise
function add_tip_pressure!(p, f, r, θ, ϕ, c, alpha, V, c0, aspect_ratio, round)

    # compute tip lift curve slope
    if aspect_ratio < 2.0
        aratio = 0.5
    elseif 2.0 <= aspect_ratio <= 24.0
        aratio = aspect_ratio_correction(aspect_ratio)
    elseif aspect_ratio > 24.0
        aratio = 1.0
    end

    # modify tip angle of attack using tip lift curve slope
    alpha = aratio*alpha

    # Mach Number
    M = V/c0

    # directivity function
    Dh = Dhfunc(M, θ, ϕ)

    if round
        l = 0.008*alpha*c
    else
        if abs(alpha) <= 2.0
            l = (0.0230 + 0.0169*alpha)*c
        else
            l = (0.0378 + 0.0095*alpha)*c
        end
    end

    # Maximum Mach Number (eq. 64)
    Mmax = (1.0 + 0.036*alpha)*M

    # Maximum velocity
    Vmax = Mmax * c0

    # compute scale
    term = M^2*Mmax^3*l^2*Dh/r^2
    if iszero(term)
        scale = 0.0
    else
        scale = 10.0 * log10(term)
    end

    # compute pressure for each frequency
    for i = 1:length(f)

        St = (f[i]*l)/Vmax

        spl_tip = 126.0 - 30.5*(log10(St) + 0.3)^2 + scale

        p[i] += 10^(spl_tip/10)
    end

    return p
end

function boundary_thickness(Re, c, alpha)

    d0 = c*10.0^(1.6569 - 0.9045*log10(Re) + 0.0596*log10(Re)^2)

    dp = d0*10^(-0.04175*alpha+0.00106*alpha^2)

    return dp
end

function displacement_thickness(Re, c, alpha, trip)

    if trip # heavily tripped boundary layer

        # displacement thickness (eq 3)
        if Re <= 0.3e6
            d0 = c*0.0601*Re^(-0.114)
        else
            d0 = c*10.0^(3.411 - 1.5397*log10(Re) + 0.1059*(log10(Re))^2)
        end
        # itrip == 2 && d0 *= 0.6

    else # untripped boundary layer

        # displacement thickness (eq 6)
        d0 = c*10.0^(3.0187 - 1.5397*log10(Re) + 0.1059*(log10(Re))^2)

    end

    # --- Pressure Side Displacement Thickness --- #

    # displacement thickness (eq 9)
    dp = d0*10.0^(-0.0432*alpha + 0.00113*alpha^2)
    # itrip == 3 && dp *= 1.48

    # --- Suction Side Displacement Thickness --- #

    if trip # heavily tripped boundary layer

        # displacement thickness on suction side (eq 12)
        if alpha <= 5.0
            ds = d0*10.0^(0.0679*alpha)
        elseif 5.0 < alpha <= 12.5
            ds = d0*0.381*10.0^(0.1516*alpha)
        else # 12.5 < alpha
            ds = d0*14.296*10.0^(0.0258*alpha)
        end

    else # untripped boundary layer

        # displacement thickness on suction side (eq 15)
        if alpha <= 7.5
            ds = d0*10.0^(0.0679*alpha)
        elseif 7.5 < alpha <= 12.5
            ds = d0*0.0162*10.0^(0.3066*alpha)
        else # 12.5 < alpha
            ds = d0*52.42*10.0^(0.0258*alpha)
        end

    end

    return dp, ds
end

end # module
