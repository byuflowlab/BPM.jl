# BPM

[![Build Status](https://travis-ci.org/moore54/BPM.jl.svg?branch=master)](https://travis-ci.org/moore54/BPM.jl)

[![Coverage Status](https://coveralls.io/repos/moore54/BPM.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/moore54/BPM.jl?branch=master)

[![codecov.io](http://codecov.io/github/moore54/BPM.jl/coverage.svg?branch=master)](http://codecov.io/github/moore54/BPM.jl?branch=master)



# BPM Turbine Acoustics

Turbine acoustic code using the BPM equations developed by Brooks, Pope, and Marcolini

Developed by Eric Tingey at Brigham Young University, 2015-2017,
Translated to Julia by Kevin Moore and Taylor McDonald at BYU, 2018

This code models the acoustic propagation of a wind turbine based on turbulent boundary layer edge noise, separation stall noise, tip vortex formation noise, laminar boundary layer vortex shedding noise, and trailing edge bluntness vortex shedding noise. Turbulent inflow noise is not assumed in this current code. The semi-empirical equations were developed from the NACA 0012 airfoil data and the blade segments used in the test file are based on the NREL 5-MW wind turbine. Scaling of the segments is based on the blade length specified.

Brooks, T., Pope, D., and Marcolini, M., “Aipower Self-Noise and Prediction,” NASA, 1989.

Brooks, T., and Marcolini, M., "Airfoil Tip Vortex Formation Noise," AIAA Journal, 1986.

Vargas, L., "Wind Turbine Noise Prediction," Master's Thesis, Technical University of Lisbon, 2008.


## Installation instructions

Pkg.clone("url-to-package-on-git")

## Running the Julia code

This julia code can be run from another file using:
```julia

import BPM

SPL_HAWT = BPM.turbinepos(turbx, turby, obs, winddir, windvel, rpm, B, h, rad,\
c, c1, alpha, nu, c0, psi, AR, noise_corr)
```
"""
Calculating the sound pressure level for a HAWT

Parameters
----------
turbx : array
    x-positions of all the turbines heard by an observer (east to west, meter)
turby : array
    y-positions of all the turbines heard by an observer (north to south, meter)
obs : array
    x-, y-, and z-position of a specified observer (E-W, N-S, height; meter)
winddir : float
    direction the wind blows from (180=N, 270=E, 0=S, 90=W; degree)
windvel : array
    wind velocity at each turbine in the specified wind direction (m/s)
rpm : array
    rotation rate of each turbine (RPM)
B : float
    number of blades on a turbine
h : float
    height of a turbine (meter)
rad : array
    radial positions of the blade geometry (meter)
c : array
    chord length at each radial segment (meter)
c1 : array
    distance from the pitch axis to leading edge at each radial segment (meter)
alpha : array
    angle of attack of each radial segment (degree)
nu : float
    kinematic viscosity of the air (m^2/s)
c0 : float
    speed of sound of the air (m/s)
psi : float
    solid angle of turbine blades between upper and lower sides of trailing edge (degree)
AR : float
    aspect ratio of turbine blades
noise_corr : float
    correction factor for SPL calculations (1=none, use if calculations differ from expected)

Returns
----------
SPL_HAWT : float
    sound pressure level calculated at observer location (dB)
"""

## Run tests

```julia
Pkg.test("BPM")
```
