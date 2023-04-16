# BPM.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/BPM.jl/dev)
![](https://github.com/byuflowlab/BPM.jl/workflows/CI/badge.svg)

**A semi-empirical code for modeling the acoustics of propellers and turbines** 

 - Developed by Eric Tingey at the FLOW Lab in Brigham Young University, 2015-2017,
 - Translated to Julia by Kevin Moore and Taylor McDonald at FLOW Lab, 2018
 - Refactored by Tyler Critchfield and Eduardo Alvarez at FLOW Lab, 2020.
 - Completely rewritten (with added verification cases) by Taylor McDonnell in 2023.

## Summary:

This code models the acoustic propagation of propellers and wind turbines using the semi-empirical BPM equations developed by Brooks, Pope, and Marcolini.  

These equations capture the following sources of noise:
 - Laminar boundary layer vortex shedding noise
 - Turbulent boundary layer edge noise
 - Separation stall noise
 - Trailing edge bluntness vortex shedding noise  
 - Tip vortex formation noise

Note that his code does not currently model turbulent inflow noise. 

## Installation:

```julia
] add https://github.com/byuflowlab/BPM.jl
```

## Usage:

This package exports the following function:

```julia
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
```

## References

Brooks, T., Pope, D., and Marcolini, M., “Aipower Self-Noise and Prediction,” NASA, 1989.

Brooks, T., and Marcolini, M., "Airfoil Tip Vortex Formation Noise," AIAA Journal, 1986.

Vargas, L., "Wind Turbine Noise Prediction," Master's Thesis, Technical University of Lisbon, 2008.
