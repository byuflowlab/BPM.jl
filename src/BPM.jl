# The MIT License (MIT)
#
# Copyright (c) 2015 Eric Tingey, Brigham Young University
#
# Modified 2018 Kevin Moore, Taylor McDonald BYU

module BPM
# Subroutine to use the BPM equations for turbine acoustics
export turbinepos,turbinepos_VAWT
# cubic spline interpolation setup (for Tip Vortex Noise)
function splineint(n, x, y, xval)
    yval = 0.0
    # assuming the values of x are in accending order
    for i = 1:n
        if (xval < x[i])
            if (i == 2)
                x1 = x[1]
                x2 = x[2]
                x3 = x[3]
                y1 = y[1]
                y2 = y[2]
                y3 = y[3]
                yval = cubspline(x1,x2,x3,y1,y2,y3,xval)
            elseif (i == n)
                x1 = x[n-2]
                x2 = x[n-1]
                x3 = x[n]
                y1 = y[n-2]
                y2 = y[n-1]
                y3 = y[n]
                yval = cubspline(x1,x2,x3,y1,y2,y3,xval)
            else
                if (xval <= (x[i]+x[i-1])/2.0)
                    x1 = x[i-2]
                    x2 = x[i-1]
                    x3 = x[i]
                    y1 = y[i-2]
                    y2 = y[i-1]
                    y3 = y[i]
                    yval = cubspline(x1,x2,x3,y1,y2,y3,xval)
                else
                    x1 = x[i-1]
                    x2 = x[i]
                    x3 = x[i+1]
                    y1 = y[i-1]
                    y2 = y[i]
                    y3 = y[i+1]
                    yval = cubspline(x1,x2,x3,y1,y2,y3,xval)
                end
            end
            break
        elseif (xval == x[i])
            yval = y[i]
            break
        end

    end
    return yval
end

function cubspline(x1, x2, x3, y1, y2, y3, xval)

    a11 = 2.0/(x2-x1)
    a12 = 1.0/(x2-x1)
    a13 = 0.0
    a21 = 1.0/(x2-x1)
    a22 = 2.0*((1.0/(x2-x1))+(1.0/(x3-x2)))
    a23 = 1.0/(x3-x2)
    a31 = 0.0
    a32 = 1.0/(x3-x2)
    a33 = 2.0/(x3-x2)
    b1 = 3.0*(y2-y1)/(x2-x1)^2
    b2 = 3.0*(((y2-y1)/(x2-x1)^2)+((y3-y2)/(x3-x2)^2))
    b3 = 3.0*(y3-y2)/(x3-x2)^2

    bot = a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32
    if (xval < x2)
        xtop = b1*a22*a33+a12*a23*b3+a13*b2*a32-a13*a22*b3-a12*b2*a33-b1*a23*a32
        ytop = a11*b2*a33+b1*a23*a31+a13*a21*b3-a13*b2*a31-b1*a21*a33-a11*a23*b3

        k1 = xtop/bot
        k2 = ytop/bot

        a = k1*(x2-x1)-(y2-y1)
        b = -k2*(x2-x1)+(y2-y1)
        t = (xval-x1)/(x2-x1)

        yval = (1.0-t)*y1+t*y2+t*(1.0-t)*(a*(1.0-t)+b*t)
    else
        ytop = a11*b2*a33+b1*a23*a31+a13*a21*b3-a13*b2*a31-b1*a21*a33-a11*a23*b3
        ztop = a11*a22*b3+a12*b2*a31+b1*a21*a32-b1*a22*a31-a12*a21*b3-a11*b2*a32

        k2 = ytop/bot
        k3 = ztop/bot

        a = k2*(x3-x2)-(y3-y2)
        b = -k3*(x3-x2)+(y3-y2)
        t = (xval-x2)/(x3-x2)

        yval = (1.0-t)*y2+t*y3+t*(1.0-t)*(a*(1.0-t)+b*t)
    end
    return yval
end

# Function to compute directivity angles and distance
# Based on work by Luis Vargas (Wind Turbine Noise Prediction)
function direct(n, xt, yt, zt, c, c1, d, Hub, beta)

    theta_e = zeros(n)
    phi_e = zeros(n)
    c2 = zeros(n)
    r = zeros(n)

    # distance from pitch-axis to trailing edge
    c2[1:n] = c[1:n]-c1[1:n]

    # Calculating observer location from hub
    xo = xt # lateral direction
    yo = yt # downstream direction
    zo = zt-Hub # height direction

    for i=1:n
        # Calculating trailing edge position from hub
        xs = sin(beta)*d[i]-cos(beta)*c2[i]
        zs = cos(beta)*d[i]+sin(beta)*c2[i]

        # Calculating observer position from trailing edge
        xe_d = xo-xs
        ze_d = zo-zs

        # Rotating observer position with repsect to beta
        theta = pi-beta
        xe = cos(theta)*xe_d+sin(theta)*ze_d
        ze = -sin(theta)*xe_d+cos(theta)*ze_d

        # Calculating observer distance and directivity angles
        r[i] = sqrt(yo^2+xe^2+ze^2)
        theta_e[i] = atan2(sqrt(yo^2+ze^2),xe)
        phi_e[i] = atan2(yo,ze)

        # Quadratic smoothing when phi_e is close to 0 or 180 degrees
        if (abs(phi_e[i]) < 5.0*pi/180.0)
            if (phi_e[i] >= 0.0)
                sign = 1
            else
                sign = -1
            end
            phi_er = abs(phi_e[i])*180.0/pi
            phi_er = 0.1*phi_er^2+2.5
            phi_e[i] = sign*phi_er*pi/180.0
        elseif (abs(phi_e[i]) > 175.0*pi/180.0)
            if (phi_e[i] >= 0.0)
                sign = 1
            else
                sign = -1
            end
            phi_er = abs(phi_e[i])*180.0/pi
            phi_er = -0.1*(phi_er-180.0)^2+177.5
            phi_e[i] = sign*phi_er*pi/180.0
        end

    end
    return r,theta_e,phi_e
end #direct

#VAWT Directivty function
function directVAWT(n, xt, yt, zt, c, c1, ht, rad, Hub, rotdir, beta)

    theta_e = zeros(n)
    phi_e = zeros(n)
    c2 = zeros(n)
    r = zeros(n)

    # distance from pitch-axis to trailing edge
    c2[1:n] = c[1:n]-c1[1:n]

    # Calculating observer location from hub
    xo = xt # lateral direction
    yo = yt # downstream direction
    zo = zt-Hub # height direction

    for i=1:n
        # Calculating trailing edge position from hub
        if rotdir >= 0.0
            xs=-rad*cos(beta)-c2[i]*sin(beta)
            ys=-rad*sin(beta)+c2[i]*cos(beta)
            zs=ht[i]
        else
            xs=-rad*cos(beta)+c2[i]*sin(beta)
            ys=-rad*sin(beta)-c2[i]*cos(beta)
            zs=ht[i]
        end

        # Calculating observer position from trailing edge
        xe_d = xo-xs
        ye_d = yo-ys
        ze = zo-zs

        # Rotating observer position with repsect to beta
        theta = rotdir*(pi/2.0)+beta
        xe = cos(theta)*xe_d+sin(theta)*ye_d
        ye = -sin(theta)*xe_d+cos(theta)*ye_d

        # Calculating observer distance and directivity angles
        r[i] = sqrt(xe^2+ye^2+ze^2)
        theta_e[i] = atan2(sqrt(ye^2+ze^2),xe)
        phi_e[i] = atan2(ye,ze)
        # Quadratic smoothing when theta_e or phi_e are close to 0, +/-180 degrees
        if (abs(theta_e[i])< (5.0*pi/180.0))
            if (theta_e[i] >= 0.0)
                sign = 1.0
            else
                sign = -1.0
            end
            theta_er=abs(theta_e[i])*180.0/pi
            theta_er=0.1*(theta_er)^2+2.5
            theta_e[i] = sign*theta_er*pi/180.0
        elseif (abs(theta_e[i])> (175.0*pi/180))
            if (theta_e[i] >=0 )
                sign =1.0
            else
                sign=-1.0
            end
            theta_er = abs(theta_e[i])*180.0/pi
            theta_er = -0.1*(theta_er-180.0)^2+177.5
            theta_e[i] = sign*theta_er*pi/180.0
        end

        if (abs(phi_e[i]) < 5.0*pi/180.0)
            if (phi_e[i] >= 0.0)
                sign = 1.0
            else
                sign = -1.0
            end
            phi_er = abs(phi_e[i])*180.0/pi
            phi_er = 0.1*phi_er^2+2.5
            phi_e[i] = sign*phi_er*pi/180.0
        elseif (abs(phi_e[i]) > 175.0*pi/180.0)
            if (phi_e[i] >= 0.0)
                sign = 1
            else
                sign = -1
            end
            phi_er = abs(phi_e[i])*180.0/pi
            phi_er = -0.1*(phi_er-180.0)^2+177.5
            phi_e[i] = sign*phi_er*pi/180.0
        end
    end
    return r,theta_e,phi_e
end #directVAWT

# Directivity function for high-frequency noise
# not for high-angle separation; becomes inaccurate for theta_e approaching 180 deg
function Dhfunc(theta_e, phi_e,M)
    conv = 0.8 # convection factor for speed
    Mc = M*conv

    Dh = (2.0*(sin(theta_e/2.0))^2*(sin(phi_e))^2)/((1.0+M*cos(theta_e))*(1.0+(M-Mc)*cos(theta_e))^2)
    return Dh
end

# Directivity function for low-frequency noise
function Dlfunc(theta_e, phi_e,M)
    Dl = ((sin(theta_e))^2*(sin(phi_e))^2)/(1.0+M*cos(theta_e))^4
    return Dl
end #Dlfunc

# Spectral Function A
function Afunc(ain, Re)
    a = abs(log10(ain))

    # Calculating Amin
    if (a < 0.204)
        Amin = sqrt(67.552-886.788*a^2)-8.219
    elseif (a >= 0.204 && a <= 0.244)
        Amin = -32.665*a+3.981
    else
        Amin = -142.795*a^3+103.656*a^2-57.757*a+6.006
    end

    # Calculating Amax
    if (a < 0.13)
        Amax = sqrt(67.552-886.788*a^2)-8.219
    elseif (a >= 0.13 && a <= 0.321)
        Amax = -15.901*a+1.098
    else
        Amax = -4.669*a^3+3.491*a^2-16.699*a+1.149
    end

    # Calculating a0
    if (Re < 9.52e4)
        a0 = 0.57
    elseif (Re >= 9.52e4 && Re <= 8.57e5)
        a0 = -9.57e-13*(Re-8.57e5)^2+1.13
    else
        a0 = 1.13
    end

    # Calculating Amin(a0)
    if (a0 < 0.204)
        Amin0 = sqrt(67.552-886.788*a0^2)-8.219
    elseif (a0 >= 0.204 && a0 <= 0.244)
        Amin0 = -32.665*a0+3.981
    else
        Amin0 = -142.795*a0^3+103.656*a0^2-57.757*a0+6.006
    end

    # Calculating Amax(a0)
    if (a0 < 0.13)
        Amax0 = sqrt(67.552-886.788*a0^2)-8.219
    elseif (a0 >= 0.13 && a0 <= 0.321)
        Amax0 = -15.901*a0+1.098
    else
        Amax0 = -4.669*a0^3+3.491*a0^2-16.699*a0+1.149
    end

    AR0 = (-20.0-Amin0)/(Amax0-Amin0)

    Aspec = Amin+AR0*(Amax-Amin)
    return Aspec
end #Afunc

# Spectral Function B
function Bfunc(bin, Re)

    b = abs(log10(bin))

    # Calculating Bmin
    if (b < 0.13)
        Bmin = sqrt(16.888-886.788*b^2)-4.109
    elseif (b >= 0.13 && b <= 0.145)
        Bmin = -83.607*b+8.138
    else
        Bmin = -817.810*b^3+355.210*b^2-135.024*b+10.619
    end

    # Calculating Bmax
    if (b < 0.10)
        Bmax = sqrt(16.888-886.788*b^2)-4.109
    elseif (b >= 0.10 && b <= 0.187)
        Bmax = -31.330*b+1.854
    else
        Bmax = -80.541*b^3+44.174*b^2-39.381*b+2.344
    end

    # Calculating b0
    if (Re < 9.52e4)
        b0 = 0.30
    elseif (Re >= 9.52e4 && Re <= 8.57e5)
        b0 = -4.48e-13*(Re-8.57e5)^2+0.56
    else
        b0 = 0.56
    end

    # Calculating Bmin(b0)
    if (b0 < 0.13)
        Bmin0 = sqrt(16.888-886.788*b0^2)-4.109
    elseif (b0 >= 0.13 && b0 <= 0.145)
        Bmin0 = -83.607*b0+8.138
    else
        Bmin0 = -817.810*b0^3+355.210*b0^2-135.024*b0+10.619
    end

    # Calculating Bmax(b0)
    if (b0 < 0.10)
        Bmax0 = sqrt(16.888-886.788*b0^2)-4.109
    elseif (b0 >= 0.10 && b0 <= 0.187)
        Bmax0 = -31.330*b0+1.854
    else
        Bmax0 = -80.541*b0^3+44.174*b0^2-39.381*b0+2.344
    end

    BR0 = (-20.0-Bmin0)/(Bmax0-Bmin0)

    Bspec =  Bmin+BR0*(Bmax-Bmin)
    return Bspec
end #Bfunc

function G1func(e)
    if (e <= 0.5974)
        G1 = 39.8*log10(e)-11.12
    elseif (e <= 0.8545 && e > 0.5974)
        G1 = 98.409*log10(e)+2.0
    elseif (e <= 1.17 && e > 0.8545)
        G1 = sqrt(2.484-506.25*(log10(e))^2)-5.076
    elseif (e <= 1.674 && e > 1.17)
        G1 = -98.409*log10(e)+2.0
    else
        G1 = -39.8*log10(e)-11.12
    end
    return G1
end #G1func

function G2func(d)
    if (d <= 0.3237)
        G2 = 77.852*log10(d)+15.328
    elseif (d <= 0.5689 && d > 0.3237)
        G2 = 65.188*log10(d)+9.125
    elseif (d <= 1.7579 && d > 0.5689)
        G2 = -114.052*(log10(d))^2
    elseif (d <= 3.0889 && d > 1.7579)
        G2 = -65.188*log10(d)+9.125
    else
        G2 = -77.852*log10(d)+15.328
    end
    return G2
end #G2func

function G3func(alpha)
    G3 = 171.04-3.03*alpha
    return G3
end #G3func

function G4func(hdav, psi)
    if (hdav <= 5.0)
        G4 = 17.5*log10(hdav)+157.5-1.114*psi
    else
        G4 = 169.7-1.114*psi
    end
    return G4
end #G4func

function G5func(hdav, psi, StSt_peak)
    # finding G5 at phi = 14 deg
    eta = log10(StSt_peak)

    if (hdav < 0.25)
        mu = 0.1221
    elseif (hdav < 0.62 && hdav >= 0.25)
        mu = -0.2175*hdav+0.1755
    elseif (hdav < 1.15 && hdav >= 0.62)
        mu = -0.0308*hdav+0.0596
    else
        mu = 0.0242
    end

    if (hdav <= 0.02)
        m = 0.0
    elseif (hdav <= 0.5 && hdav > 0.02)
        m = 68.724*(hdav)-1.35
    elseif (hdav <= 0.62 && hdav > 0.5)
        m = 308.475*hdav-121.23
    elseif (hdav <= 1.15 && hdav > 0.62)
        m = 224.811*hdav-69.35
    elseif (hdav <= 1.2 && hdav > 1.15)
        m = 1583.28*hdav-1631.59
    else
        m = 268.344
    end

    eta_0 = -sqrt((m^2*mu^4)/(6.25+m^2*mu^2))
    k = 2.5*sqrt(1.0-(eta_0/mu)^2)-2.5-m*eta_0

    if (eta < eta_0)
        G14 = m*eta+k
    elseif (eta < 0.0 && eta >= eta_0)
        G14 = 2.5*sqrt(1.0-(eta/mu)^2)-2.5
    elseif (eta < 0.03616 && eta >= 0.0)
        G14 = sqrt(1.5625-1194.99*eta^2)-1.25
    else
        G14 = -155.543*eta+4.375
    end

    # finding G5 at psi = 0 deg
    hdav_prime = 6.724*hdav^2-4.019*hdav+1.107

    if (hdav_prime < 0.25)
        mu0 = 0.1221
    elseif (hdav_prime < 0.62 && hdav_prime >= 0.25)
        mu0 = -0.2175*hdav_prime+0.1755
    elseif (hdav_prime < 1.15 && hdav_prime >= 0.62)
        mu0 = -0.0308*hdav_prime+0.0596
    else
        mu0 = 0.0242
    end

    if (hdav_prime <= 0.02)
        m0 = 0.0
    elseif (hdav_prime <= 0.5 && hdav_prime > 0.02)
        m0 = 68.724*hdav_prime-1.35
    elseif (hdav_prime <= 0.62 && hdav_prime > 0.5)
        m0 = 308.475*hdav_prime-121.23
    elseif (hdav_prime <= 1.15 && hdav_prime > 0.62)
        m0 = 224.811*hdav_prime-69.35
    elseif (hdav_prime <= 1.2 && hdav_prime > 1.15)
        m0 = 1583.28*hdav_prime-1631.59
    else
        m0 = 268.344
    end

    eta_00 = -sqrt((m0^2*mu0^4)/(6.25+m0^2*mu0^2))
    k0 = 2.5*sqrt(1.0-(eta_00/mu0)^2)-2.5-m0*eta_00

    if (eta < eta_00)
        G0 = m0*eta+k0
    elseif (eta < 0.0 && eta >= eta_00)
        G0 = 2.5*sqrt(1.0-(eta/mu0)^2)-2.5
    elseif (eta < 0.03616 && eta >= 0.0)
        G0 = sqrt(1.5625-1194.99*eta^2)-1.25
    else
        G0 = -155.543*eta+4.375
    end

    G5 = G0+0.0714*psi*(G14-G0)
    return G5
end #G5func

# Turbulent Boundary Layer Trailing Edge Noise
function TBLTEfunc(f, V, L, c, r, theta_e, phi_e, alpha, nu, c0, trip)
    # constants
    M = V/c0
    Re = (V*c)/nu

    if (trip == false)
        # UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
        d0 = c*(10.0^(1.6569-0.9045*log10(Re)+0.0596*(log10(Re))^2))
        d0_d = c*(10.0^(3.0187-1.5397*log10(Re)+0.1059*(log10(Re))^2))
    else
        # TRIPPED boundary layer at 0 deg- thickness, displacement thickness
        d0 = c*(10.0^(1.892-0.9045*log10(Re)+0.0596*(log10(Re))^2))
        if (Re <= 0.3e6)
            d0_d = c*0.0601*Re^(-0.114)
        else
            d0_d = c*(10.0^(3.411-1.5397*log10(Re)+0.1059*(log10(Re))^2))
        end
    end

    # boundary layer on pressure side- thickness, displacement thickness
    dpr = d0*(10.0^(-0.04175*alpha+0.00106*alpha^2))
    dp_d = d0_d*(10.0^(-0.0432*alpha+0.00113*alpha^2))

    if (trip == false)
        # UNTRIPPED boundary layer on suction side- displacement thickness
        if (alpha <= 7.5 && alpha >= 0.0)
            ds_d = d0_d*10.0^(0.0679*alpha)
        elseif (alpha <= 12.5 && alpha > 7.5)
            ds_d = d0_d*0.0162*10.0^(0.3066*alpha)
            #elseif (alpha <= 25.0 && alpha > 12.5)
        else
            ds_d = d0_d*52.42*10.0^(0.0258*alpha)
        end
    else
        # TRIPPED boundary layer on suction side- displacement thickness
        if (alpha <= 5.0 && alpha >= 0.0)
            ds_d = d0_d*10.0^(0.0679*alpha)
        elseif (alpha <= 12.5 && alpha > 5.0)
            ds_d = d0_d*0.381*10.0^(0.1516*alpha)
            #elseif (alpha <= 25.0 && alpha > 12.5)
        else
            ds_d = d0_d*14.296*10.0^(0.0258*alpha)
        end
    end

    Dh = Dhfunc(theta_e,phi_e,M)
    Dl = Dlfunc(theta_e,phi_e,M)

    Stp = (f*dp_d)/V
    Sts = (f*ds_d)/V

    St1 = 0.02*M^(-0.6)

    if (alpha < 1.33)
        St2 = St1*1.0
    elseif (alpha <= 12.5 && alpha >= 1.33)
        St2 = St1*10.0^(0.0054*(alpha-1.33)^2)
    else
        St2 = St1*4.72
    end

    St_bar = (St1+St2)/2.0

    St_peak = max(St1,St2,St_bar)

    apre = Stp/St1
    asuc = Sts/St1
    bang = Sts/St2

    gamma = 27.094*M+3.31
    gamma0 = 23.43*M+4.651
    beta = 72.65*M+10.74
    beta0 = -34.19*M-13.82

    if (Re < 2.47e5)
        K1 = -4.31*log10(Re)+156.3
    elseif (Re >= 2.47e5 && Re <= 8.0e5)
        K1 = -9.0*log10(Re)+181.6
    else
        K1 = 128.5
    end

    if (alpha < (gamma0-gamma))
        K2 = K1-1000.0
    elseif (alpha >= (gamma0-gamma) && alpha <= (gamma0+gamma))
        K2 = K1+sqrt(beta^2-(beta/gamma)^2*(alpha-gamma0)^2)+beta0
    else
        K2 = K1-12.0
    end

    Re_dp = (V*dp_d)/nu

    if (Re_dp <= 5000.0)
        DeltaK1 = alpha*(1.43*log10(Re_dp)-5.29)
    else
        DeltaK1 = 0.0
    end

    # Keeping observer distance from getting too close to the turbine
    if (r < 1e-8)
        rc = 1e-8
    else
        rc = r
    end

    if (alpha > 12.5 || alpha > gamma0)
        # Turbulent Boundary Layer Separation Stall Noise (TBLSS); this is were the airfoil is stalling and stall noise dominates
        # SPLp = -infinity; 10^(SPLp/10) = 0
        # SPLs = -infinity; 10^(SPLs/10) = 0

        A = Afunc(bang,3.0*Re)

        SPLa = 10.0*log10((ds_d*M^5*L*Dl)/rc^2)+A+K2

        TBLTE = 10.0*log10(10.0^(SPLa/10.0))

    else
        Ap = Afunc(apre,Re)
        As = Afunc(asuc,Re)
        B = Bfunc(bang,Re)

        SPLp = 10.0*log10((dp_d*M^5*L*Dh)/rc^2)+Ap+(K1-3.0)+DeltaK1
        SPLs = 10.0*log10((ds_d*M^5*L*Dh)/rc^2)+As+(K1-3.0)
        SPLa = 10.0*log10((ds_d*M^5*L*Dh)/rc^2)+B+K2

        TBLTE =  10.0*log10(10.0^(SPLp/10.0)+10.0^(SPLs/10.0)+10.0^(SPLa/10.0))

    end
    return TBLTE
end #TBLTEfunc

# Turbulent Boundary Layer Tip Vortex Noise
function TBLTVfunc(f, V, c, r, theta_e, phi_e, atip, c0, tipflat, AR)

    # constants
    M = V/c0
    Mmax = M*(1.0+0.036*atip)

    Dh = Dhfunc(theta_e,phi_e,M)

    # Tip vortex noise correction based on data from "Airfoil Tip Vortex Formation Noise"
    AR_data = [2.0,2.67,4.0,6.0,12.0,24.0]
    atipcorr_data = [0.54,0.62,0.71,0.79,0.89,0.95]

    if ((AR >= 2.0) && (AR <= 24.0))
        atipcorr = splineint(6,AR_data,atipcorr_data,AR)
    elseif (AR > 24.0)
        atipcorr = 1.0
    else
        atipcorr = 0.5
    end

    atip_d = atip*atipcorr

    if (tipflat == false)
        # rounded tip
        l = 0.008*c*atip_d
    else
        # flat tip
        if (atip_d <= 2.0 && atip_d >= 0.0)
            l = c*(0.0230+0.0169*atip_d)
        else
            l = c*(0.0378+0.0095*atip_d)
        end
    end

    St = (f*l)/(V*(1.0+0.036*atip_d))

    # Keeping observer distance from getting too close to the turbine
    if (r < 1e-8)
        rc = 1e-8
    else
        rc = r
    end

    TBLTV = 10.0*log10((M^2*Mmax^3*l^2*Dh)/rc^2)-30.5*(log10(St)+0.3)^2+126.0
    return TBLTV
end #TBLTVfunc

# Laminar Boundary Layer Vortex Shedding
function LBLVSfunc(f, V, L, c, r, theta_e, phi_e, alpha, nu, c0, trip)

    # constants
    M = V/c0
    Re = (V*c)/nu

    if (trip == false)
        # UNTRIPPED boundary layer at 0 deg- thickness
        d0 = c*(10.0^(1.6569-0.9045*log10(Re)+0.0596*(log10(Re))^2))
    else
        # TRIPPED boundary layer at 0 deg- thickness
        d0 = c*(10.0^(1.892-0.9045*log10(Re)+0.0596*(log10(Re))^2))
    end
    # boundary layer on pressure side- thickness
    dpr = d0*(10.0^(-0.04175*alpha+0.00106*alpha^2))

    St = (f*dpr)/V

    Dh = Dhfunc(theta_e,phi_e,M)

    if (Re <= 1.3e5)
        St1 = 0.18
    elseif (Re <= 4.0e5 && Re > 1.3e5)
        St1 = 0.001756*Re^0.3931
    else
        St1 = 0.28
    end

    St_peak = St1*10.0^(-0.04*alpha)

    e = St/St_peak

    G1 = G1func(e)

    if (alpha <= 3.0)
        Re0 = 10.0^(0.215*alpha+4.978)
    else
        Re0 = 10.0^(0.12*alpha+5.263)
    end

    d = Re/Re0

    G2 = G2func(d)
    G3 = G3func(alpha)

    # Keeping observer distance from getting too close to the turbine
    if (r < 1e-8)
        rc = 1e-8
    else
        rc = r
    end

    LBLVS = 10.0*log10((dpr*M^5*L*Dh)/rc^2)+G1+G2+G3
    return LBLVS
end #LBLVSfunc

# Trailing Edge Bluntness Vortex Shedding Noise
function TEBVSfunc(f, V, L, c, h, r, theta_e, phi_e, alpha, nu, c0, psi, trip)

    # constants
    M = V/c0
    Re = (V*c)/nu

    if (trip == false)
        # UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
        d0 = c*(10.0^(1.6569-0.9045*log10(Re)+0.0596*(log10(Re))^2))
        d0_d = c*(10.0^(3.0187-1.5397*log10(Re)+0.1059*(log10(Re))^2))
    else
        # TRIPPED boundary layer at 0 deg- thickness, displacement thickness
        d0 = c*(10.0^(1.892-0.9045*log10(Re)+0.0596*(log10(Re))^2))
        if (Re <= 0.3e6)
            d0_d = c*0.0601*Re^(-0.114)
        else
            d0_d = c*(10.0^(3.411-1.5397*log10(Re)+0.1059*(log10(Re))^2))
        end
    end

    # boundary layer on pressure side- thickness, displacement thickness
    dpr = d0*(10.0^(-0.04175*alpha+0.00106*alpha^2))
    dp_d = d0_d*(10.0^(-0.0432*alpha+0.00113*alpha^2))

    if (trip == false)
        # UNTRIPPED boundary layer on suction side- displacement thickness
        if (alpha <= 7.5 && alpha >= 0.0)
            ds_d = d0_d*10.0^(0.0679*alpha)
        elseif (alpha <= 12.5 && alpha > 7.5)
            ds_d = d0_d*0.0162*10.0^(0.3066*alpha)
            #elseif (alpha <= 25 && alpha > 12.5)
        else
            ds_d = d0_d* 52.42* 10.0^(0.0258*alpha)
        end
    else
        # TRIPPED boundary layer on suction side- displacement thickness
        if (alpha <= 5.0 && alpha >= 0.0)
            ds_d = d0_d*10.0^(0.0679*alpha)
        elseif (alpha <= 12.5 && alpha > 5.0)
            ds_d = d0_d*0.381*10.0^(0.1516*alpha)
            #elseif (alpha <= 25.0 && alpha > 12.5)
        else
            ds_d = d0_d*14.296*10.0^(0.0258*alpha)
        end
    end

    Dh = Dhfunc(theta_e,phi_e,M)
    St = (f*h)/V
    dav = (dp_d+ds_d)/2.0

    hdav = h/dav

    if (hdav >= 0.2)
        St_peak = (0.212-0.0045*psi)/(1.0+0.235*(1.0/hdav)-0.0132*(1.0/hdav)^2)
    else
        St_peak = 0.1*(hdav)+0.095-0.00243*psi
    end

    StSt_peak = St/St_peak

    G4 = G4func(hdav,psi)
    G5 = G5func(hdav,psi,StSt_peak)

    # Keeping observer distance from getting too close to the turbine
    if (r < 1e-8)
        rc = 1e-8
    else
        rc = r
    end

    TEBVS = 10.0*log10((h*M^(5.5)*L*Dh)/rc^2)+G4+G5
    return TEBVS
end #TEBVSfunc

# Computing the overall sound pressure level (OASPL) of a turbine defined below (in dB)
function OASPL(ox, oy, oz, windvel, rpm, B, Hub, rad, c, c1, alpha, nu, c0, psi, AR)
    # constants

    nf = 27
    bf = 3

    n = length(rad)

    L = zeros(n-1)
    d = zeros(n-1)
    V = zeros(n-1)
    h = zeros(n-1)
    TV_t = zeros(B)
    TE_t = zeros((n-1)*B)
    BLVS_t = zeros((n-1)*B)
    BVS_t = zeros((n-1)*B)

    TE = zeros(nf)
    TV = zeros(nf)
    BLVS = zeros(nf)
    BVS = zeros(nf)
    SPLf = zeros(bf, nf)
    SPLfA = zeros(bf, nf)
    SPLoa_d = zeros(bf)


    # Using untripped or tripped boundary layer specficiation
    trip = false # untripped
    # trip = true # tripped

    # Tip specfication
    tipflat = false # round
    # tipflat = true # flat

    # Parameters of the wind turbine
    omega = (rpm*2.0*pi)/60.0  # angular velocity (rad/sec)

    for i = 1:n-1
        L[i] = rad[i+1]-rad[i] # length of each radial section (m)
        d[i] = rad[i] # radial section to be used in directivity calculations (m)
        V[i] = sqrt((omega*rad[i])^2+windvel^2) # wind speed over the blade (m/s)
    end

    h[1:n-1] = 0.01*c[1:n-1]  # trailing edge thickness; 1% of chord length (m)
    atip = alpha[n-1]  # angle of attack of the tip region (deg)

    # Blade rotation increments to rotate around (45 deg from Vargas paper)
    # beta = [0.0,0.25*pi,0.5*pi,0.75*pi,pi,1.25*pi,1.5*pi,1.75*pi] # 8 increments
    beta = [0.0,2.0*pi/9.0,4.0*pi/9.0] # 3 increments (equivalent of 9 for 3 blades)
    # beta = [0.0,pi] # 2 increments
    # beta = [0.0] # 1 increment (top blade facing straight up)

    B_int = 2.0*pi/B # Intervals between blades (from the first blade at 0 deg)

    # One-third octave band frequencies (Hz)
    f = [100.0,125.0,160.0,200.0,250.0,315.0,400.0,500.0,
    630.0,800.0,1000.0,1250.0,1600.0,2000.0,2500.0,3150.0,
    4000.0,5000.0,6300.0,8000.0,10000.0,12500.0,16000.0,
    20000.0,25000.0,31500.0,40000.0]

    # A-weighting curve (dBA) for sound perception correction
    AdB = [-19.145,-16.190,-13.244,-10.847,-8.675,-6.644,
    -4.774,-3.248,-1.908,-0.795,0.0,0.576,0.993,1.202,
    1.271,1.202,0.964,0.556,-0.114,-1.144,-2.488,-4.250,
    -6.701,-9.341,-12.322,-15.694,-19.402]

    for di=1:bf # for each rotation increment
        for j=1:nf # for each frequency
            for bi=1:B # for each blade
                # Calcuating observer distances and directivty angles for the given blade orientation
                r,theta_e,phi_e = direct(n-1,ox,oy,oz,c,c1,d,Hub,beta[di]+(bi-1)*B_int)

                TBLTV = TBLTVfunc(f[j],V[n-1],c[n-1],r[n-1],theta_e[n-1],phi_e[n-1],atip,c0,
                tipflat,AR)
                TV_t[bi] = TBLTV
                for k=1:n-1
                    # Calculating sound pressure level (dB) for each noise source at each radial position
                    TBLTE = TBLTEfunc(f[j],V[k],L[k],c[k],r[k],theta_e[k],phi_e[k],alpha[k],
                    nu,c0,trip)

                    if (trip == false)
                        LBLVS = LBLVSfunc(f[j],V[k],L[k],c[k],r[k],theta_e[k],phi_e[k],alpha[k],
                        nu,c0,trip)
                    else
                        LBLVS = 0.0
                    end
                    TEBVS = TEBVSfunc(f[j],V[k],L[k],c[k],h[k],r[k],theta_e[k],phi_e[k],
                    alpha[k],nu,c0,psi,trip)

                    # Assigning noise to blade segment
                    TE_t[k+[n-1]*(bi-1)] = TBLTE
                    BLVS_t[k+[n-1]*(bi-1)] = LBLVS
                    BVS_t[k+[n-1]*(bi-1)] = TEBVS
                end
            end

            # Adding sound pressure levels (dB)
            TE[j] = 10.0*log10(sum(10.0.^(TE_t/10.0)))
            TV[j] = 10.0*log10(sum(10.0.^(TV_t/10.0)))
            BLVS[j] = 10.0*log10(sum(10.0.^(BLVS_t/10.0)))
            BVS[j] = 10.0*log10(sum(10.0.^(BVS_t/10.0)))

            # Combining noise sources into overall SPL
            SPLf[di,j] = 10.0*log10(10.0^(TE[j]/10.0)+10.0^(TV[j]/10.0)+10.0^(BLVS[j]/10.0)+10.0^(BVS[j]/10.0))
        end


        # Correcting with A-weighting
        # SPLf[1:nf] = SPLf[1:nf]+AdB[1:nf]
        SPLfA[di,1:nf] = SPLf[di,1:nf]+AdB[1:nf]

        # Adding SPLs for each rotation increment
        SPLoa_d[di] = 10.0*log10(sum(10.0.^(SPLfA[di,:]/10.0)))

        # Protecting total calcuation from negative SPL values
        if (SPLoa_d[di] < 0.0)
            SPLoa_d[di] = 0.0
        end
    end

    # Performing root mean square calculation of SPLs at rotation increments for final value
    SPLoa = sqrt(sum(SPLoa_d.^2)/bf)
    SPLf_overall = sqrt.(sum(SPLf.^2, 1)/bf)
    SPLfA_overall = sqrt.(sum(SPLfA.^2, 1)/bf)

    return SPLoa, SPLf_overall, SPLfA_overall
end #OASPL

#Computing the overall sound pressure level (OASPL) of a turbine defined below (in dB)
function OASPLVAWT(p, ox, oy, oz, B, Hub, high, rad, c, c1, alpha, nu, 
                    c0, psi, rot, Vinf, wakex, wakey, AR)

   # constants
   nf = 27
   bf = 8

   n = length(high)
   L = zeros(n-1)
   d = zeros(n-1)
   V = zeros(n-1)
   h = zeros(n-1)
   TV_t = zeros(2*B)
   TE_t = zeros((n-1)*B)
   BLVS_t = zeros((n-1)*B)
   BVS_t = zeros((n-1)*B)

   TE = zeros(nf)
   TV = zeros(nf)
   BLVS = zeros(nf)
   BVS = zeros(nf)
   SPLf = zeros(nf)
   SPLoa_d = zeros(bf)
   theta_vel = zeros(p)
   highmid = zeros(n-1)

   # Using untripped or tripped boundary layer specficiation
   trip = false # untripped
   # trip = true # tripped

   # Tip specfication
   tipflat = false # round
   # tipflat = true # flat
   for i = 1:n-1
       L[i] = high[i+1]-high[i] # length of each height section (m)
       highmid[i] = (high[i+1]+high[i])/2.0
   end
   h[1:n-1] = 0.01*c[1:n-1]
   atip1 = alpha[1] #angle of attack of the tip region on bottom (deg)
   atip2= alpha[n-1] #angle of attack of the tip region on top (deg)
   if (rot >= 0.0)
       rotdir = 1.0
   else
       rotdir = -1.0
   end

   # Blade rotation increments to rotate around (45 deg from Vargas paper)
   beta = [0.0,0.25*pi,0.5*pi,0.75*pi,pi,1.25*pi,1.5*pi,1.75*pi] # 8 increments
   #beta = [0.0,2.0*pi/9.0,4.0*pi/9.0] # 3 increments (equivalent of 9 for 3 blades)
   # beta = [0.0,pi] # 2 increments
   # beta = [0.0] # 1 increment (top blade facing straight up)

   B_int = 2.0*pi/B # Intervals betwqeen blades (from the first blade at 0 deg)

   # One-third octave band frequencies (Hz)
   f = [100.0,125.0,160.0,200.0,250.0,315.0,400.0,500.0,
   630.0,800.0,1000.0,1250.0,1600.0,2000.0,2500.0,3150.0,
   4000.0,5000.0,6300.0,8000.0,10000.0,12500.0,16000.0,
   20000.0,25000.0,31500.0,40000.0]

   # A-weighting curve (dBA) for sound perception correction
   AdB = [-19.145,-16.190,-13.244,-10.847,-8.675,-6.644,
   -4.774,-3.248,-1.908,-0.795,0.0,0.576,0.993,1.202,
   1.271,1.202,0.964,0.556,-0.114,-1.144,-2.488,-4.250,
   -6.701,-9.341,-12.322,-15.694,-19.402]

   for i=1:p
       theta_vel[i] = (2.0*pi/p)*i-(2.0*pi/p)/2.0
   end

   for di=1:bf # for each rotation increment
       for j=1:nf # for each frequency
           for bi=1:B # for each blade
               # Calcuating observer distances and directivty angles for the given blade orientation
               theta = beta[di]+(bi-1)*B_int
               r,theta_e,phi_e = directVAWT(n-1,ox,oy,oz,c,c1,highmid,rad,Hub,rotdir,theta)
               if ((theta >= theta_vel[1]) & (theta <= theta_vel[p]))
                   velwx=splineint(p,theta_vel,wakex,theta)
                   velwy=splineint(p,theta_vel,wakey,theta)
               else
                   velwx = (wakex[1]+wakex[p])/2.0
                   velwy = (wakey[1]+wakey[p])/2.0
               end

               Vx=rot*rad*cos(theta)+Vinf+velwx
               Vy=rot*rad*sin(theta) + velwy
               V=sqrt(Vx^2+Vy^2)
               TBLTV = TBLTVfunc(f[j],V[1],c[1],r[1],theta_e[1],phi_e[1],atip1,c0,
               tipflat,AR)
               TV_t[2*(bi-1)+1] = TBLTV
               TBLTV= TBLTVfunc(f[j],V,c[n-1],r[n-1],theta_e[n-1],phi_e[n-1],atip2,c0,
               tipflat,AR)
               TV_t[2*(bi-1)+2] = TBLTV

               for k=1:n-1
                   # Calculating sound pressure level (dB) for each noise source at each radial position
                   TBLTE = TBLTEfunc(f[j],V,L[k],c[k],r[k],theta_e[k],phi_e[k],alpha[k],
                   nu,c0,trip)

                   if (trip == false)
                       LBLVS = LBLVSfunc(f[j],V,L[k],c[k],r[k],theta_e[k],phi_e[k],alpha[k],
                       nu,c0,trip)
                   else
                       LBLVS = 0.0
                   end
                   TEBVS = TEBVSfunc(f[j],V,L[k],c[k],h[k],r[k],theta_e[k],phi_e[k],
                   alpha[k],nu,c0,psi,trip)
                   # Assigning noise to blade segment
                   TE_t[k+[n-1]*(bi-1)] = TBLTE
                   BLVS_t[k+[n-1]*(bi-1)] = LBLVS
                   BVS_t[k+[n-1]*(bi-1)] = TEBVS
               end
           end
           # Adding sound pressure levels (dB)
           TE[j] = 10.0*log10(sum(10.0.^(TE_t/10.0)))
           TV[j] = 10.0*log10(sum(10.0.^(TV_t/10.0)))
           BLVS[j] = 10.0*log10(sum(10.0.^(BLVS_t/10.0)))
           BVS[j] = 10.0*log10(sum(10.0.^(BVS_t/10.0)))

           # Combining noise sources into overall SPL
           SPLf[j] = 10.0*log10(10.0^(TE[j]/10.0)+10.0^(TV[j]/10.0)+10.0^(BLVS[j]/10.0)+10.0^(BVS[j]/10.0))
       end


       # Correcting with A-weighting
       SPLf[1:nf] = SPLf[1:nf]+AdB[1:nf]

       # Adding SPLs for each rotation increment
       SPLoa_d[di] = 10.0*log10(sum(10.0.^(SPLf/10.0)))

       # Protecting total calcuation from negative SPL values
       if (SPLoa_d[di] < 0.0)
           SPLoa_d[di] = 0.0
       end
   end

   # Performing root mean square calculation of SPLs at rotation increments for final value
   SPLoa = sqrt(sum(SPLoa_d.^2)/bf)
   return SPLoa
end #OASPLVAWT

# Placing a turbine in a specified location and finding the OASPL of the turbine with reference to an observer
"""

turbinepos(x,y,obs,winddir,windvel,rpm,B,Hub,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr)

Calculating the sound pressure level for a HAWT

# Parameters
----------
- `turbx::array`:  x-positions of all the turbines heard by an observer (east to west, meter)
- `turby::array`:  y-positions of all the turbines heard by an observer (north to south, meter)
- `obs::array`:  x-, y-, and z-position of a specified observer (E-W, N-S, height; meter)
- `winddir::float`:  direction the wind blows from (180=N, 270=E, 0=S, 90=W; degree)
- `windvel::array`:  wind velocity at each turbine in the specified wind direction (m/s)
- `rpm::array`:  rotation rate of each turbine (RPM)
- `B::float`:  number of blades on a turbine
- `h::float`:  height of a turbine (meter)
- `rad::array`:  radial positions of the blade geometry (meter)
- `c::array`:  chord length at each radial segment (meter)
- `c1::array`:  distance from the pitch axis to leading edge at each radial segment (meter)
- `alpha::array`:  angle of attack of each radial segment (degree)
- `nu::float`:  kinematic viscosity of the air (m^2/s)
- `c0::float`:  speed of sound of the air (m/s)
- `psi::float`:  solid angle of turbine blades between upper and lower sides of trailing edge (degree)
- `AR::float`:  aspect ratio of turbine blades
- `noise_corr::float`:  correction factor for SPL calculations (1=none, use if calculations differ from expected)

# Returns
----------
- `SPL_HAWT::float`:  sound pressure level calculated at observer location (dB)
"""
function turbinepos(x, y, obs, winddir, windvel, rpm, 
                    B, Hub, rad, c, c1, alpha, nu, c0, 
                    psi, AR, noise_corr)

    nturb = length(x)
    tSPL = zeros(nturb)
    SPLf = zeros(nturb, 27) #based on nf from OASPL function
    SPLfA = zeros(nturb, 27)
    windrad = (winddir+180.0)*pi/180.0

    for i = 1:nturb # for each turbine
        # Centering the turbine at (0,0) with repect to the observer location
        ox = obs[1]-x[i]
        oy = obs[2]-y[i]
        oz = obs[3]

        # Adjusting the coordinates to turbine reference frame (wind moving in y-direction)
        rxy = sqrt(ox^2+oy^2)
        ang = atan2(oy,ox)+windrad

        ox = rxy*cos(ang)
        oy = rxy*sin(ang)

        # Calculating the overall SPL of each of the turbines at the observer location
        tSPL[i], SPLf[i,:], SPLfA[i,:] = OASPL(ox,oy,oz,windvel[i],rpm[i],B,Hub,rad,c,c1,alpha,nu,c0,psi,AR)
    end

    # Combining the SPLs from each turbine and correcting the value based on the wind farm
    SPL_obs = (10.0*log10(sum(10.0.^(tSPL/10.0)))) * noise_corr
    SPLf = (10.0*log10(sum(10.0.^(SPLf/10.0), 1))) * noise_corr
    SPLfA = (10.0*log10(sum(10.0.^(SPLfA/10.0), 1))) * noise_corr

    return SPL_obs, SPLf, SPLfA
end #turbinepos
"""

turbinepos_VAWT(p,x,y,obs,winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,Vinf,wakex,wakey)

Calculating the sound pressure level for a HAWT

# Parameters
----------
- `p::integer`: number of points along blade flight path to calculate velocities
- `turbx::array`:  x-positions of all the turbines heard by an observer (east to west, meter)
- `turby::array`:  y-positions of all the turbines heard by an observer (north to south, meter)
- `obs::array`:  x-, y-, and z-position of a specified observer (E-W, N-S, height; meter)
- `winddir::float`:  direction the wind blows from (180=N, 270=E, 0=S, 90=W; degree)
- `B::float`:  number of blades on a turbine
- `Hub::float`: hub height of a turbine (meter)
- `high::array`: height positions along the turbine blade (meter)
- `rad::array`:  turbine radius (meter)
- `c::array`:  chord length at each radial segment (meter)
- `c1::array`:  distance from the pitch axis to leading edge at each radial segment (meter)
- `alpha::array`:  angle of attack of each radial segment (degree)
- `nu::float`:  kinematic viscosity of the air (m^2/s)
- `c0::float`:  speed of sound of the air (m/s)
- `psi::float`:  solid angle of turbine blades between upper and lower sides of trailing edge (degree)
- `AR::float`:  aspect ratio of turbine blades
- `noise_corr::float`:  correction factor for SPL calculations (1=none, use if calculations differ from expected)
- `rot::array`: rotation rate of each turbine (rad/s)
- `Vinf::float`: free stream wind speed (m/s)
- `wakex::array`: the wake influenced x-velcoity of the turbine at each point along the blade flight path (m/s)
- `wakey::array`: the wake influenced y-velcoity of the turbine at each point along the blade flight path (m/s)
# Returns
----------
- `SPL_VAWT::float`:  sound pressure level calculated at observer location (dB)
"""
function turbinepos_VAWT(p, x, y, obs, winddir, B, Hub, high,
                            rad, c, c1, alpha, nu, c0, psi, AR, 
                            noise_corr, rot, Vinf, wakex, wakey)

    nturb = length(x)
    tSPL = zeros(nturb)
    # SPLf = zeros(nturb)
    # SPLfA = zeros(nturb)
    
    wakexd = zeros(p)
    wakeyd = zeros(p)
    windrad = (winddir+180.0)*pi/180.0
    for i = 1:nturb # for each turbine
        # Centering the turbine at (0,0) with repect to the observer location
        ox = obs[1]-x[i]
        oy = obs[2]-y[i]
        oz = obs[3]

        # Adjusting the coordinates to turbine reference frame (wind moving in y-direction)
        rxy = sqrt(ox^2+oy^2)
        ang = atan2(oy,ox)+windrad

        ox = rxy*cos(ang)
        oy = rxy*sin(ang)

        for j=1:p
            k=p*(i-1)+j
            wakexd[j] = wakex[k]
            wakeyd[j] = wakey[k]
        end
        # Calculating the overall SPL of each of the turbines at the observer location
        tSPL[i] = OASPLVAWT(p,ox,oy,oz,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,rot[i],Vinf,wakexd,wakeyd,AR)
    end

    # Combining the SPLs from each turbine and correcting the value based on the wind farm
    SPL_obs = (10.0*log10(sum(10.0.^(tSPL/10.0))))*noise_corr

    return SPL_obs
end #turbinepos_VAWT

end # module
