# Project Reorganization Changes

This document outlines the key changes made to reorganize the Monte Carlo simulations codebase for CrI2 and CrI3 monolayers.

## Directory Structure Changes

1. **Modular Source Code Organization**:
   - Organized Fortran source files into logical modules:
     - `src/core/`: Core simulation components
     - `src/io/`: Input/output operations
     - `src/simulation/`: Simulation drivers
     - `src/tools/`: Utility tools

2. **Separate Data and Examples**:
   - Created `data/` directory with subdirectories for different simulation types
   - Added `examples/` directory with example input files for each material

3. **Improved Build System**:
   - Implemented a unified Makefile for all simulation components
   - Standardized compilation flags and build targets

4. **Script Organization**:
   - Moved Python scripts to dedicated `scripts/` directory
   - Improved script interfaces with proper command-line arguments

5. **Documentation**:
   - Added `docs/` directory with detailed project documentation

## Code Improvements

1. **Modernized Run Script**:
   - Replaced interactive command-line queries with command-line arguments
   - Added proper parameter validation
   - Implemented better error handling
   - Added support for reading parameters from input files

2. **Enhanced Data Analysis**:
   - Improved data loading and processing
   - Added plotting capabilities
   - Added support for combining and analyzing data from multiple runs

3. **Better Code Structure**:
   - Clearly separated interface between Fortran simulation code and Python analysis tools
   - Standardized input/output file formats

## Documentation Additions

1. **Project README**:
   - Comprehensive overview of the project
   - Installation and usage instructions
   - Description of key parameters and components

2. **Code Documentation**:
   - Detailed description of project structure
   - Module and function descriptions
   - Examples of how to run simulations

## Benefits of the New Structure

1. **Easier Maintenance**:
   - Clear separation of concerns
   - Modular components with well-defined interfaces
   - Consistent file organization

2. **Better Extensibility**:
   - Can easily add new material types
   - Can implement new analysis methods
   - Can extend with additional simulation types

3. **Improved Usability**:
   - Clear examples for getting started
   - Consistent command-line interfaces
   - Better organized output data

4. **Reproducibility**:
   - Standard input format ensures reproducible simulations
   - Clear connection between input parameters and output data 