# Sarin Heavy Gas Dispersion Simulation

This repository contains a Julia implementation of a lightweight Lagrangian particle model for rapid assessment of Sarin vapor dispersion in an urban environment.

The `sarin.jl` script simulates the movement of heavy toxic gas particles while accounting for:

- gravity-driven settling of dense vapor;
- turbulent particle dispersion;
- basement entrapment and stagnation zones;
- ventilation-driven particle escape from basements;
- stack-effect-driven vertical infiltration in buildings;
- Sarin mass degradation due to hydrolysis and solar photodegradation;
- time-dependent evaporation from a liquid spill source.

The model is intended for fast scenario-based computational experiments and decision-support prototyping. It does not solve the full Navier--Stokes equations, but instead uses simplified physical parameterizations suitable for near real-time hazard assessment.

## Requirements

- Julia 1.9 or newer
- Standard Julia libraries only:
  - `Random`
  - `Printf`

No external packages are required.

## Running the Simulation

Run the script from the command line:

```bash
julia sarin.jl
```

The program prints the calculated degradation constants and time-dependent estimates of active airborne mass, mass trapped in basements, and the number of simulated particles.

## Main File

The full source code is available here:

https://github.com/iamvadimov/gasdyn/blob/main/sarin.jl

## Notes

This code is a research-oriented demonstration model. Physical constants, kinetic parameters, source characteristics, and geometry settings should be reviewed and adjusted before applying the model to any real-world scenario.
