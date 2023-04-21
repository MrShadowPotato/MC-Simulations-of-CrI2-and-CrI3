import pandas as pd
import numpy as np


class SimulationData:
    directory = 'remote/stabilizer/'

    def __init__(self, n, t, mcs, orientation, seeds=list) -> None:
        self.n = n
        self.t = t
        self.mcs = mcs
        self.orientation = orientation
        self.seeds = seeds
        self.data_frames = []
        self.energy_arrays = []
        self.magnetization_arrays = []
        self.data_frames_sum = []


    def set_seeds(self, seeds):
        if len(self.seeds) == 0:
            self.seeds = seeds
        else: 
            print('Seeds already set')
            print('Seeds: ', self.seeds)
    
    def get_data(self, nlines_to_read):
        #Will read the data from the files and append the data to the data_frames list.
        #Will also append the energy and magnetization arrays to the energy_arrays and magnetization_arrays lists.
        ncommentary_lines = 3
        start_index = self.mcs - nlines_to_read + ncommentary_lines
        for seed in self.seeds:
            path = f'{self.directory}n{self.n}t{self.t}m{self.mcs}s{seed}o{self.orientation}.txt'
            df = pd.read_table(path, header=None, names=['step','energy','magnetization'], delim_whitespace=True, skiprows=start_index, skipfooter=1, engine='python')
            self.data_frames.append(df)

            array = np.genfromtxt(path, skip_header=3, skip_footer=1, usecols=(0,1,2), unpack=True)
            #step_array = array[0]
            energy_array = array[1]
            magnetization_array = array[2]
            
            self.energy_arrays.append(energy_array)
            self.magnetization_arrays.append(magnetization_array)

    def write_seed_sum(self, number_of_seeds):
        #Will write the sum of the first number_of_seeds to a file.
        #And will also append the data to the data_frames_sum list.
        combined_energy = np.sum(self.energy_arrays[:number_of_seeds], axis=0)/number_of_seeds
        combined_magnetization = np.sum(self.magnetization_arrays[:number_of_seeds], axis=0)/number_of_seeds
        step_array = np.arange(1, len(combined_energy) + 1)
        
        
        file_name = f'{self.directory}n{self.n}t{self.t}m{self.mcs}number_of_seeds{number_of_seeds}o{self.orientation}.txt'
        with open(file_name, 'w') as file:
            file.write(f'# n = {self.n}, t = {self.t}, mcs = {self.mcs}, orientation = {self.orientation}, number of seeds = {number_of_seeds}\n')
            file.write('# step energy magnetization\n')
            np.savetxt(file, np.c_[step_array, combined_energy, combined_magnetization], fmt='%d %f %f')

seeds = [-111, -222, -333, -555]#-444, -555]
print('Seeds: ', seeds)

n = 80 #input('n:' )
t = 5 #input('t:' ) 
mcs = 5000 #input('mcs:' )
orientation = 2 #input('orientation:' )
start_index = int(input('Start index:' ))


instance = SimulationData(n, t, mcs, orientation, seeds)
instance.get_data(1500)
for i in range(1, len(seeds) + 1):
    instance.write_seed_sum(i)
for count, frame in enumerate(instance.data_frames_sum):
    print(f'Values for {count + 1} seeds')
    print('Energy mean', frame['energy'].mean())
    print('Energy sem', frame['energy'].sem())

        
