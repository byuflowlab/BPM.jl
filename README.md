# BPM Turbine Acoustics

Turbine acoustic code using the BPM equations developed by Brooks, Pope, and Marcolini

Developed by Eric Tingey at the FLOW Lab in Brigham Young University, 2015-2017,
Translated to Julia by Kevin Moore and Taylor McDonald at FLOW Lab, 2018
Refactored by Tyler Critchfield and Eduardo Alvarez at FLOW Lab, 2020.

This code models the acoustic propagation of a wind turbine based on turbulent boundary layer edge noise, separation stall noise, tip vortex formation noise, laminar boundary layer vortex shedding noise, and trailing edge bluntness vortex shedding noise. Turbulent inflow noise is not assumed in this current code. The semi-empirical equations were developed from the NACA 0012 airfoil data and the blade segments used in the test file are based on the NREL 5-MW wind turbine. Scaling of the segments is based on the blade length specified.

Brooks, T., Pope, D., and Marcolini, M., “Aipower Self-Noise and Prediction,” NASA, 1989.

Brooks, T., and Marcolini, M., "Airfoil Tip Vortex Formation Noise," AIAA Journal, 1986.

Vargas, L., "Wind Turbine Noise Prediction," Master's Thesis, Technical University of Lisbon, 2008.


## Running the code

This julia code can be run from another file using:
```julia
import BPM

OASPL_HAWT, SPLf_HAWT, SPLfA_HAWT = BPM.turbinepos(turbx, turby, obs, winddir, windvel, rpm, B, h, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)

OASPL_VAWT = BPM.turbinepos_VAWT(p,x,y,obs,winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,Vinf,wakex,wakey)
```

## Input and Output Definition
```julia
"""

turbinepos(x,y,obs,winddir,windvel,rpm,B,Hub,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr)


Calculating the sound pressure level for a HAWT

Parameters
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

Returns
----------
- `OASPL_HAWT::float`:  A-weighted overall sound pressure level calculated at observer location (dB)
- `SPLf_HAWT`: sound pressure level at each frequency
- `SPLfA_HAWT`: A-weighted sound pressure level at each frequency

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
```
