# Monte Carlo Simulations of Magnetic Properties in CrI2 and CrI3 Monolayers

This repository contains Fortran-based Monte Carlo simulation code for studying the magnetic properties of two-dimensional CrI2 and CrI3 monolayers. The simulations use the Metropolis algorithm to explore the magnetic phase transitions, anisotropy effects, and field-dependent behaviors in these materials.

## Overview

2D chromium iodide materials (CrI2 and CrI3) exhibit interesting magnetic properties that can be modeled using Monte Carlo simulations. This code implements:

- Heisenberg exchange interactions
- Single-ion magnetic anisotropy
- External magnetic field effects (Zeeman interaction)
- Temperature and field-dependent simulations

## Project Structure

```
.
├── src/                    # Source code
│   ├── core/               # Core simulation components
│   │   ├── constants.f90   # Physical constants and system parameters
│   │   ├── rng.f90         # Random number generation
│   │   ├── energy.f90      # Energy calculations
│   │   ├── metropolis.f90  # Metropolis algorithm implementation
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
│   ├── input.CrI2          # Example input for CrI2
│   ├── input.CrI3          # Example input for CrI3
├── data/                   # Data output directory
│   ├── tLoop/              # Temperature dependence results
│   ├── hLoop/              # Field dependence results
├── docs/                   # Documentation
```

## Installation

### Prerequisites

- Fortran compiler (gfortran or ifort)
- Python 3.x with numpy, pandas, and matplotlib

### Building

```bash
cd src
make all
```

## Usage

1. Configure the simulation parameters in an input file (see examples directory)
2. Run a simulation:

```bash
# Using Python runner (recommended)
python scripts/run.py --program tLoop --input examples/input.CrI3

# Direct execution
cd src
./tLoop < ../examples/input.CrI3
```

3. Analyze the results:

```bash
python scripts/data_analysis.py --data-dir data/tLoop --output plots
```

## Simulation Parameters

The simulation takes various parameters including:
- System size (nx, ny): Dimensions of the 2D lattice
- Compound: Material type (CrI2 or CrI3)
- J: Exchange coupling strength
- K: Anisotropy constant
- H: External magnetic field
- T: Temperature
- Monte Carlo steps (mcs): Number of simulation steps
- Equilibration steps (neq): Steps for equilibration

## Key Files

- `energy.f90`: Calculates the energy of the system, including exchange interactions, anisotropy, and Zeeman effects
- `metropolis.f90`: Implements the Metropolis algorithm for spin updates
- `tLoop.f90`: Performs temperature-dependent simulations
- `hLoop.f90`: Performs magnetic field-dependent simulations

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this code for academic work, please cite:

```
@article{YourName2023,
  title={Monte Carlo Simulations of Magnetic Properties in CrI2 and CrI3 Monolayers},
  author={Your Name},
  journal={Your University/Institution},
  year={2023}
}
```




