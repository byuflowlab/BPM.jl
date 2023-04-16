# BPM.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/BPM.jl/dev)
![](https://github.com/byuflowlab/BPM.jl/workflows/CI/badge.svg)

**A semi-empirical acoustics code for modeling the acoustic propagation of propellers and turbines**

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

```@docs
sound_pressure_levels
```

## References

Brooks, T., Pope, D., and Marcolini, M., “Aipower Self-Noise and Prediction,” NASA, 1989.

Brooks, T., and Marcolini, M., "Airfoil Tip Vortex Formation Noise," AIAA Journal, 1986.

Vargas, L., "Wind Turbine Noise Prediction," Master's Thesis, Technical University of Lisbon, 2008.
