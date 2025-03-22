import os
import sys
import subprocess
import argparse
from pathlib import Path

# Define default parameters
DEFAULT_PARAMS = {
    'compound': 'CrI3',  # CrI3 or CrI2
    'nx': 10,            # Number of unit cells in x direction
    'ny': 10,            # Number of unit cells in y direction
    'J': 2.5,            # Exchange constant
    'K': 0.0,            # Anisotropy constant
    'iH': 0.0,           # Initial magnetic field
    'iT': 1.0,           # Initial temperature
    'g': 2.0,            # g-factor
    'mcs': 5000,         # Monte Carlo steps
    'neq': 1000,         # Equilibration steps
    'seed': 12345,       # Random seed
    'spins_orientation': 1,  # 1: All random, 2: Ordered
    'mag_dir': 'easy',   # Initial magnetization direction
    'H_vector_direction': 'easy',  # Magnetic field direction
    'fT': 50.0,          # Final temperature (for tLoop)
    'dT': 1.0,           # Temperature step (for tLoop)
    'dH': 0.5,           # Magnetic field step (for hLoop)
    'dS': 0.1,           # Spin multiplier for new spin proposals
}

# Seeds for multiple runs
seeds = [12345, 67890, 13579, 24680, 11223]

def output_file(program, iseed, n):
    """Generate output filename based on parameters"""
    return f"{program}_{DEFAULT_PARAMS['compound']}_n{n}_seed{iseed}"

def write_input(program, params, output_file_name):
    """Write input parameters to the input file"""
    with open('input', 'w') as f:
        f.write('Name of the file where the output data will be stored:\n')
        f.write(f'{output_file_name}\n')
        f.write('**** System parameters ****\n')
        f.write(f'nx= {params["nx"]}\n')
        f.write(f'ny= {params["ny"]}\n')
        f.write(f'compound= {params["compound"]}\n')
        f.write(f'J= {params["J"]}\n')
        f.write(f'K= {params["K"]}\n')
        f.write(f'iH= {params["iH"]}\n')
        f.write(f'iT= {params["iT"]}\n')
        f.write(f'g= {params["g"]}\n')
        f.write('**** Montecarlo parameters ****\n')
        f.write(f'mcs= {params["mcs"]}\n')
        f.write(f'neq= {params["neq"]}\n')
        f.write(f'seed= {params["seed"]}\n')
        f.write('**** Initial spin config ****\n')
        f.write(f'order= {params["spins_orientation"]}\n')
        f.write(f'mag_dir= {params["mag_dir"]}\n')
        f.write('**** H_vector_direction ****\n')
        f.write(f'h_dir= {params["H_vector_direction"]}\n')
        f.write('******************************************\n')
        f.write('Data for tLoop \n')
        f.write(f'fT= {params["fT"]}\n')
        f.write(f'dT= {params["dT"]}\n')
        f.write('Data for hLoop \n')
        f.write(f'dH= {params["dH"]}\n')
        f.write('**** Spin multiplier ****\n')
        f.write(f'dS= {params["dS"]}')

def parse_input_file(input_file):
    """Parse parameters from an input file"""
    params = DEFAULT_PARAMS.copy()
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found.")
        sys.exit(1)
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
        
    # Skip the first line (output filename)
    i = 1
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('nx='):
            params['nx'] = int(line.split('=')[1].strip())
        elif line.startswith('ny='):
            params['ny'] = int(line.split('=')[1].strip())
        elif line.startswith('compound='):
            params['compound'] = line.split('=')[1].strip()
        elif line.startswith('J='):
            params['J'] = float(line.split('=')[1].strip())
        elif line.startswith('K='):
            params['K'] = float(line.split('=')[1].strip())
        elif line.startswith('iH='):
            params['iH'] = float(line.split('=')[1].strip())
        elif line.startswith('iT='):
            params['iT'] = float(line.split('=')[1].strip())
        elif line.startswith('g='):
            params['g'] = float(line.split('=')[1].strip())
        elif line.startswith('mcs='):
            params['mcs'] = int(line.split('=')[1].strip())
        elif line.startswith('neq='):
            params['neq'] = int(line.split('=')[1].strip())
        elif line.startswith('seed='):
            params['seed'] = int(line.split('=')[1].strip())
        elif line.startswith('order='):
            params['spins_orientation'] = int(line.split('=')[1].strip())
        elif line.startswith('mag_dir='):
            params['mag_dir'] = line.split('=')[1].strip()
        elif line.startswith('h_dir='):
            params['H_vector_direction'] = line.split('=')[1].strip()
        elif line.startswith('fT='):
            params['fT'] = float(line.split('=')[1].strip())
        elif line.startswith('dT='):
            params['dT'] = float(line.split('=')[1].strip())
        elif line.startswith('dH='):
            params['dH'] = float(line.split('=')[1].strip())
        elif line.startswith('dS='):
            params['dS'] = float(line.split('=')[1].strip())
        i += 1
    
    return params

def run_simulation(program, params, seed_index=0, skip_compile=False):
    """Run the simulation with the given parameters"""
    # Set the seed for this run
    params['seed'] = seeds[seed_index] if seed_index < len(seeds) else params['seed']
    
    # Define paths
    project_root = Path(__file__).parent.parent
    bin_dir = project_root / "bin"
    program_path = bin_dir / program
    
    # Check if executable exists
    if not program_path.exists():
        if skip_compile:
            print(f"Error: Executable {program_path} not found. Please compile first.")
            sys.exit(1)
        else:
            print(f"Executable {program_path} not found. Compiling...")
            compile_simulation()
    
    # Generate output filename
    output_name = output_file(program, params['seed'], params['nx'])
    
    # Write input file
    write_input(program, params, output_name)
    
    print(f"Running {program} simulation with seed {params['seed']}...")
    print(f"Output will be saved as: {output_name}")
    
    # Run the simulation
    try:
        subprocess.run([str(program_path)], check=True)
        print(f"Simulation completed successfully!")
        return True
    except subprocess.CalledProcessError:
        print(f"Error running simulation.")
        return False

def compile_simulation():
    """Compile the simulation programs"""
    project_root = Path(__file__).parent.parent
    src_dir = project_root / "src"
    
    try:
        os.chdir(src_dir)
        print("Compiling simulation programs...")
        subprocess.run(["make", "all"], check=True)
        print("Compilation successful!")
        os.chdir(project_root)
        return True
    except subprocess.CalledProcessError:
        print("Error compiling simulation programs.")
        os.chdir(project_root)
        return False

def main():
    """Main entry point for the simulation runner"""
    parser = argparse.ArgumentParser(description='Run Monte Carlo simulations for CrI2/CrI3 monolayers')
    parser.add_argument('--program', choices=['tLoop', 'hLoop', 'sLoop', 'generator'], 
                       required=True, help='Simulation program to run')
    parser.add_argument('--input', help='Input file with parameters')
    parser.add_argument('--seeds', type=int, default=1, help='Number of seeds to run (max 5)')
    parser.add_argument('--compile', action='store_true', help='Compile before running')
    parser.add_argument('--skip-compile', action='store_true', help='Skip compilation check')
    
    args = parser.parse_args()
    
    # Compile if requested
    if args.compile:
        if not compile_simulation():
            sys.exit(1)
    
    # Get parameters
    params = DEFAULT_PARAMS.copy()
    if args.input:
        params = parse_input_file(args.input)
    
    # Run simulations with different seeds
    num_seeds = min(args.seeds, len(seeds))
    for i in range(num_seeds):
        if not run_simulation(args.program, params, i, args.skip_compile):
            print(f"Failed on seed {i+1}/{num_seeds}")
            break
        print(f"Completed seed {i+1}/{num_seeds}")

if __name__ == "__main__":
    main()
        