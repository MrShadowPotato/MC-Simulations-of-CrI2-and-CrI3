import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from pathlib import Path


class SimulationDataAnalyzer:
    """Class for analyzing and processing Monte Carlo simulation data"""
    
    def __init__(self, data_dir="../data", output_dir="../output"):
        """Initialize the data analyzer with directories"""
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        
        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
        self.data_frames = []
        self.energy_arrays = []
        self.magnetization_arrays = []
        self.data_frames_sum = []
    
    def set_seeds(self, seeds):
        """Set the seeds to use for analysis"""
        self.seeds = seeds
    
    def get_data(self, file_pattern, nlines_to_read=-1, skiprows=3):
        """
        Read data from files matching the pattern.
        
        Args:
            file_pattern: Pattern for files to read
            nlines_to_read: Number of lines to read (-1 for all)
            skiprows: Number of header rows to skip
        """
        files = list(self.data_dir.glob(file_pattern))
        
        if not files:
            print(f"No files found matching pattern: {file_pattern}")
            return
            
        print(f"Found {len(files)} files matching pattern {file_pattern}")
        
        for filepath in files:
            try:
                # For reading the last nlines_to_read lines
                if nlines_to_read > 0:
                    with open(filepath, 'r') as f:
                        lines = f.readlines()
                    
                    # Skip header and read only last nlines_to_read
                    start_index = max(skiprows, len(lines) - nlines_to_read)
                    
                    # Create a DataFrame from the lines
                    lines_to_read = lines[start_index:-1] if lines[-1].strip() == '' else lines[start_index:]
                    data = [line.strip().split() for line in lines_to_read]
                    
                    df = pd.DataFrame(data, columns=['step', 'energy', 'magnetization'], dtype=float)
                    self.data_frames.append(df)
                    
                    # Also store as numpy arrays
                    array = np.genfromtxt(filepath, skip_header=skiprows, skip_footer=1, usecols=(0, 1, 2), unpack=True)
                    
                    # If nlines_to_read is specified, only take last nlines_to_read
                    if nlines_to_read > 0 and len(array[0]) > nlines_to_read:
                        energy_array = array[1][-nlines_to_read:]
                        magnetization_array = array[2][-nlines_to_read:]
                    else:
                        energy_array = array[1]
                        magnetization_array = array[2]
                else:
                    # Read the entire file
                    df = pd.read_table(filepath, header=None, names=['step', 'energy', 'magnetization'], 
                                       delim_whitespace=True, skiprows=skiprows, skipfooter=1, engine='python')
                    self.data_frames.append(df)
                    
                    array = np.genfromtxt(filepath, skip_header=skiprows, skip_footer=1, usecols=(0, 1, 2), unpack=True)
                    energy_array = array[1]
                    magnetization_array = array[2]
                
                self.energy_arrays.append(energy_array)
                self.magnetization_arrays.append(magnetization_array)
                
                print(f"Successfully read data from {filepath}")
                
            except Exception as e:
                print(f"Error reading file {filepath}: {e}")
    
    def combine_data(self, number_of_files=None, output_file=None):
        """
        Combine data from multiple files into one averaged dataset
        
        Args:
            number_of_files: Number of files to combine (None for all)
            output_file: Name of output file (None for no file output)
        """
        if not self.energy_arrays:
            print("No data to combine. Run get_data() first.")
            return
        
        if number_of_files is None:
            number_of_files = len(self.energy_arrays)
        else:
            number_of_files = min(number_of_files, len(self.energy_arrays))
        
        # Calculate averages
        combined_energy = np.sum(self.energy_arrays[:number_of_files], axis=0) / number_of_files
        combined_magnetization = np.sum(self.magnetization_arrays[:number_of_files], axis=0) / number_of_files
        step_array = np.arange(1, len(combined_energy) + 1)
        
        # Create a DataFrame
        combined_df = pd.DataFrame({
            'step': step_array,
            'energy': combined_energy,
            'magnetization': combined_magnetization
        })
        
        self.data_frames_sum.append(combined_df)
        
        # Write to file if requested
        if output_file:
            output_path = self.output_dir / output_file
            with open(output_path, 'w') as file:
                file.write(f"# Combined data from {number_of_files} files\n")
                file.write("# step energy magnetization\n")
                np.savetxt(file, np.c_[step_array, combined_energy, combined_magnetization], fmt='%d %f %f')
            print(f"Combined data saved to {output_path}")
        
        return combined_df
    
    def plot_energy(self, output_file=None, title="Energy vs Steps"):
        """Plot energy vs steps for all data frames"""
        if not self.data_frames and not self.data_frames_sum:
            print("No data to plot. Run get_data() or combine_data() first.")
            return
        
        plt.figure(figsize=(10, 6))
        
        if self.data_frames:
            for i, df in enumerate(self.data_frames):
                plt.plot(df['step'], df['energy'], label=f"File {i+1}")
        
        if self.data_frames_sum:
            for i, df in enumerate(self.data_frames_sum):
                plt.plot(df['step'], df['energy'], linewidth=2, label=f"Combined {i+1}")
        
        plt.xlabel('Steps')
        plt.ylabel('Energy')
        plt.title(title)
        plt.legend()
        plt.grid(True)
        
        if output_file:
            plt.savefig(self.output_dir / output_file, dpi=300)
            print(f"Energy plot saved to {self.output_dir / output_file}")
        
        plt.show()
    
    def plot_magnetization(self, output_file=None, title="Magnetization vs Steps"):
        """Plot magnetization vs steps for all data frames"""
        if not self.data_frames and not self.data_frames_sum:
            print("No data to plot. Run get_data() or combine_data() first.")
            return
        
        plt.figure(figsize=(10, 6))
        
        if self.data_frames:
            for i, df in enumerate(self.data_frames):
                plt.plot(df['step'], df['magnetization'], label=f"File {i+1}")
        
        if self.data_frames_sum:
            for i, df in enumerate(self.data_frames_sum):
                plt.plot(df['step'], df['magnetization'], linewidth=2, label=f"Combined {i+1}")
        
        plt.xlabel('Steps')
        plt.ylabel('Magnetization')
        plt.title(title)
        plt.legend()
        plt.grid(True)
        
        if output_file:
            plt.savefig(self.output_dir / output_file, dpi=300)
            print(f"Magnetization plot saved to {self.output_dir / output_file}")
        
        plt.show()
    
    def get_statistics(self):
        """Calculate mean and standard error for energy and magnetization"""
        stats = []
        
        if self.data_frames:
            for i, df in enumerate(self.data_frames):
                stats.append({
                    'File': i+1,
                    'Energy Mean': df['energy'].mean(),
                    'Energy SEM': df['energy'].sem(),
                    'Magnetization Mean': df['magnetization'].mean(),
                    'Magnetization SEM': df['magnetization'].sem(),
                })
        
        if self.data_frames_sum:
            for i, df in enumerate(self.data_frames_sum):
                stats.append({
                    'File': f'Combined {i+1}',
                    'Energy Mean': df['energy'].mean(),
                    'Energy SEM': df['energy'].sem(),
                    'Magnetization Mean': df['magnetization'].mean(),
                    'Magnetization SEM': df['magnetization'].sem(),
                })
        
        return pd.DataFrame(stats)


def main():
    """Main entry point for data analysis"""
    parser = argparse.ArgumentParser(description='Analyze Monte Carlo simulation data')
    parser.add_argument('--data-dir', default='../data', help='Directory containing data files')
    parser.add_argument('--output-dir', default='../output', help='Directory for output files')
    parser.add_argument('--pattern', default='*.txt', help='Pattern for data files to analyze')
    parser.add_argument('--combine', type=int, default=None, help='Number of files to combine (None for all)')
    parser.add_argument('--lines', type=int, default=-1, help='Number of lines to read from end of file (-1 for all)')
    parser.add_argument('--output', default=None, help='Base name for output files')
    
    args = parser.parse_args()
    
    analyzer = SimulationDataAnalyzer(args.data_dir, args.output_dir)
    analyzer.get_data(args.pattern, args.lines)
    
    if args.combine is not None:
        output_file = f"{args.output}_combined.txt" if args.output else None
        analyzer.combine_data(args.combine, output_file)
    
    # Generate plots
    if args.output:
        analyzer.plot_energy(f"{args.output}_energy.png")
        analyzer.plot_magnetization(f"{args.output}_magnetization.png")
    else:
        analyzer.plot_energy()
        analyzer.plot_magnetization()
    
    # Print statistics
    stats = analyzer.get_statistics()
    print("\nStatistics:")
    print(stats)
    
    if args.output:
        stats.to_csv(f"{args.output_dir}/{args.output}_stats.csv", index=False)
        print(f"Statistics saved to {args.output_dir}/{args.output}_stats.csv")


if __name__ == "__main__":
    main()

        
