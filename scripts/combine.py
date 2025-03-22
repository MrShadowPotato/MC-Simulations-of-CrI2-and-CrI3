import os
import pandas as pd
import argparse
from pathlib import Path

def combine_data_files(directory, output_dir, pattern, output_file=None):
    """
    Combine data files matching the pattern into a single averaged dataset.
    
    Args:
        directory: Directory containing the data files
        output_dir: Directory to save the output
        pattern: File pattern to match
        output_file: Output file name (optional)
    
    Returns:
        Combined DataFrame if successful, None otherwise
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Get all files matching the pattern
    data_dir = Path(directory)
    files = list(data_dir.glob(pattern))
    
    if not files:
        print(f"No files found matching pattern: {pattern}")
        return None
    
    print(f"Found {len(files)} files matching {pattern}")
    
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()
    processed_files = []
    
    # Iterate over the files
    for filepath in files:
        print(f"Reading data from '{filepath}'...")
        
        try:
            # Read the data from the file
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Check if file has a proper format
            data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
            
            if not data_lines:
                print(f"Warning: No data found in {filepath}, skipping")
                continue
            
            # Check if the data has even number of rows (expected for our format)
            if len(data_lines) % 2 != 0:
                # Standard format data file
                df = pd.read_table(filepath, delim_whitespace=True, comment='#')
            else:
                # No header file
                df = pd.read_table(filepath, header=None, delim_whitespace=True, comment='#')
            
            if combined_data.empty:
                combined_data = df
            else:
                # Check if dataframes have the same shape
                if df.shape != combined_data.shape:
                    print(f"Warning: {filepath} has different structure, skipping")
                    continue
                # Add the data
                combined_data = combined_data.add(df)
            
            processed_files.append(filepath)
            
        except Exception as e:
            print(f"Error processing file {filepath}: {e}")
    
    # Calculate the average
    if processed_files:
        combined_data = combined_data.div(len(processed_files))
        print(f"Successfully combined {len(processed_files)} files.")
        
        # Save the combined data if output_file is provided
        if output_file:
            output_path = Path(output_dir) / output_file
            combined_data.to_csv(output_path, sep=' ', index=False, header=False)
            print(f"Combined data saved to '{output_path}'.")
        
        return combined_data
    else:
        print(f"No valid files were processed.")
        return None


def main():
    """Main entry point for the script"""
    parser = argparse.ArgumentParser(description='Combine Monte Carlo simulation data files')
    parser.add_argument('--data-dir', default='../data', help='Directory containing data files')
    parser.add_argument('--output-dir', default='../output', help='Directory for output files')
    parser.add_argument('--pattern', required=True, help='File pattern to match (e.g., "tLoop_CrI3_n*")')
    parser.add_argument('--output', help='Output file name')
    
    args = parser.parse_args()
    
    if not args.output:
        # Generate default output name based on pattern
        args.output = f"combined_{args.pattern.replace('*', 'X')}.txt"
    
    # Combine the data
    combined_data = combine_data_files(args.data_dir, args.output_dir, args.pattern, args.output)
    
    if combined_data is not None:
        # Print some statistics
        print("\nStatistics for combined data:")
        for col in combined_data.columns:
            print(f"{col}: Mean = {combined_data[col].mean():.6f}, Std = {combined_data[col].std():.6f}")


if __name__ == "__main__":
    main()

# Example usage

for n in range(20, 121, 10):
    print(f"Combining data for n = {n}...")
    combine_data_files('./remote/tLoop', 'tCrI3', n, 'dT0.5o2mdeasyk0H0.0mcs50000')
#for n in range(100, 121, 20):
#    print(f"Combining data for n = {n}...")
#    combine_data_files('../remote/tLoop', 'tCrI2', n, 'dT0.5o2mdeasyk0H0.0mcs50000', 5)

#sCrI3n60t1o2k0H0s5
#tCrI2n30dT1.0o2mdeasyk0H0.0mcs5000s1

#tCrI2n80dT0.5o2mdeasyk0H0.0mcs50000s1

#tCrI2n30CrI2n30lesstime?mcs5000s1