# Project Structure Documentation

## Overview

This project implements Monte Carlo simulations for studying magnetic properties of 2D chromium iodide monolayers (CrI2 and CrI3). The code is organized into modules that handle various aspects of the simulation, from energy calculations to I/O operations.

## Directory Structure

```
.
├── src/                    # Source code
│   ├── core/               # Core simulation components
│   │   ├── constants.f90   # Physical constants and system parameters
│   │   ├── rng.f90         # Random number generation
│   │   ├── energy.f90      # Energy calculations
│   │   ├── metropolis.f90  # Metropolis algorithm implementation
│   │   ├── mag.f90         # Magnetization calculations
│   │   └── ...             # Additional core modules
│   ├── io/                 # Input/output routines
│   │   ├── read_input.f90  # Input parameter handling
│   ├── simulation/         # Simulation drivers
│   │   ├── tLoop.f90       # Temperature loop simulation
│   │   ├── hLoop.f90       # Magnetic field loop simulation
│   │   ├── sLoop.f90       # Structure loop simulation
│   ├── tools/              # Utility tools
│   │   ├── generator.f90   # Structure generator for CrI2 and CrI3
│   ├── Makefile            # Unified build system
├── scripts/                # Analysis scripts
│   ├── run.py              # Simulation runner
│   ├── combine.py          # Data combination utilities
│   ├── data_analysis.py    # Analysis and plotting tools
├── examples/               # Example input files and configurations
│   ├── CrI2/               # Examples for CrI2
│   ├── CrI3/               # Examples for CrI3
├── data/                   # Data output directory
│   ├── tLoop/              # Temperature dependence results
│   ├── hLoop/              # Field dependence results
│   ├── sLoop/              # Structure dependence results
├── docs/                   # Documentation
```

## Module Descriptions

### Core Modules

1. **constants.f90**: Defines physical constants and system parameters used throughout the simulation.
   - Physical constants: Boltzmann constant, Bohr magneton, etc.
   - Global variables: system energy, spins, neighbors, etc.

2. **rng.f90**: Random number generation for the Monte Carlo simulation.
   - Random number generator for spin updates
   - Random vector generation

3. **energy.f90**: Energy calculation routines.
   - Heisenberg exchange interactions
   - Magnetic anisotropy
   - Zeeman energy (external field)

4. **metropolis.f90**: Implementation of the Metropolis algorithm.
   - Spin update proposals
   - Acceptance/rejection logic
   - System energy and magnetization updates

5. **mag.f90**: Magnetization calculations and related functions.
   - Magnetization vector calculations
   - Sublattice magnetization for antiferromagnetic systems

### I/O Module

- **read_input.f90**: Handles reading and parsing simulation parameters from input files.
  - System parameters (nx, ny, compound, etc.)
  - Monte Carlo parameters (steps, equilibration, etc.)
  - Initial conditions and field parameters

### Simulation Drivers

1. **tLoop.f90**: Temperature loop simulation.
   - Performs simulation across a range of temperatures
   - Calculates thermodynamic properties: specific heat, susceptibility

2. **hLoop.f90**: Magnetic field loop simulation.
   - Performs simulations across a range of external field values
   - Calculates field-dependent magnetization and energy

3. **sLoop.f90**: Structure loop simulation.
   - Performs simulations for different system sizes or structures
   - Useful for finite-size scaling analysis

### Tools

- **generator.f90**: Generates initial structures for CrI2 and CrI3 monolayers.
  - Creates lattice configurations
  - Sets up neighbor relationships

## Build System

The unified Makefile in the src directory handles compilation of all components. Key targets:

- `make all`: Build all executables
- `make tLoop`: Build only the temperature loop executable
- `make hLoop`: Build only the field loop executable
- `make sLoop`: Build only the structure loop executable
- `make generator`: Build only the structure generator
- `make clean`: Clean build artifacts

## Running Simulations

Simulations can be run using the `scripts/run.py` script which:

1. Reads parameters from command line or input files
2. Prepares input files for the Fortran executables
3. Runs the appropriate simulation program
4. Manages output files

Example usage:
```bash
python scripts/run.py --program tLoop --input examples/CrI3/input.CrI3
```

## Data Analysis

Two main scripts are provided for data analysis:

1. **combine.py**: Combines data from multiple simulation runs (useful for averaging).
2. **data_analysis.py**: Analyzes simulation data and generates plots.

Example usage:
```bash
python scripts/data_analysis.py --data-dir data/tLoop --pattern "tLoop_CrI3_n10*" --output analysis_CrI3_n10
``` 