import os
import pandas as pd

def combine_data_files(directory, compound, n, output_file_name, num_seeds):
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Iterate over the seeds
    for iseed in range(num_seeds):
        # Define the filename based on the seed
        filename = f'{compound}n{n}{output_file_name}s{iseed + 1}'
        print(f"Reading data from '{filename}'...")
        
        # Check if the file exists
        if os.path.exists(os.path.join(directory, filename)):
            # Read the data from the file
            path = os.path.join(directory, filename)
            with open(path, 'r') as f:
                lines = f.readlines()
                if len(lines) % 2 != 0 :    
                    df = pd.read_table(path, delim_whitespace=True)
                else:
                    df = pd.read_table(path, header=None, delim_whitespace=True)
            if combined_data.empty:
                combined_data = df
            else: 
                combined_data = combined_data.add(df)                


    if not combined_data.empty:
        combined_data = combined_data.div(num_seeds)
    
    
    # Transpose the combined data to have each data point as a row
    
    # Save the combined data to a new file
    combined_filename = f'{compound}n{n}{output_file_name}sp'
    if not combined_data.empty:
        combined_data.to_csv(os.path.join(directory, combined_filename), header=None, sep=' ', index=False)
        print(f"Combined data saved to '{combined_filename}'.")
    else:
        print(f"No data to combine for '{combined_filename}'.")
    





# Example usage

for n in range(20, 121, 10):
    print(f"Combining data for n = {n}...")
    combine_data_files('./remote/tLoop', 'tCrI3', n, 'dT0.5o2mdeasyk0H0.0mcs50000', 5)
#for n in range(100, 121, 20):
#    print(f"Combining data for n = {n}...")
#    combine_data_files('../remote/tLoop', 'tCrI2', n, 'dT0.5o2mdeasyk0H0.0mcs50000', 5)

#sCrI3n60t1o2k0H0s5
#tCrI2n30dT1.0o2mdeasyk0H0.0mcs5000s1

#tCrI2n80dT0.5o2mdeasyk0H0.0mcs50000s1

#tCrI2n30CrI2n30lesstime?mcs5000s1