import subprocess
import time

class Data:
    def __init__(self, seed, mcs, temperature, nx, ny, spins_orientation):
        self.seed = seed
        self.mcs = mcs
        self.temperature = temperature
        self.nx = nx
        self.ny = ny
        self.spins_orientation = spins_orientation
        self.output_file = self.output_file()

    def output_file(self):
        return f'nxy{self.nx}t{self.temperature}m{self.mcs}s{self.seed}o{self.spins_orientation}.txt'
    
    def write_params(self, compuond='CrI3'):
        with open ('params.txt', 'w') as f:
            f.write(f'seed = {self.seed}\n')
            f.write(f'mcs = {self.mcs}\n')
            f.write(f'temperature = {self.temperature}\n')
            f.write(f'nx = {self.nx}\n')
            f.write(f'ny = {self.ny}\n')
            f.write(f'spins_orientation = {self.spins_orientation}\n')
            f.write(f'output_file = {self.output_file}\n')
            f.write(f'Cr_xyz = n{self.nx}Cr_{compuond}.xyz\n')
            f.write(f'CrI3_xyz = n{self.nx}Cr_{compuond}neighbors.txt\n')



seeds = [-111, -222, -333, -444, -555]

choice =input('Enter 1 to run xsimulation or 2 to run xgenerate: ')

#Inputs
initial_nxy = int(input('Enter the initial nx and ny (nx = ny): '))
final_nxy = int(input('Enter the final nx and ny (nx = ny): '))
mcs = int(input('Enter the number of Monte Carlo steps: '))
temperature = int(input('Enter the temperature: '))
print('Enter the spins orientation: ')
spins_orientation = int(input('1 for all random and 2 for all the same:'))


#Create nohup
'''
'''
if choice == '1':
    for nxy in range(initial_nxy, final_nxy + 1):
        for seed in seeds:
            Data(seed, mcs, temperature, nxy, nxy, spins_orientation).write_params()
            subprocess.run(f'nohup ./xsimulation > ./stabilizer_output/n{nxy}t{temperature}s{seed}.output &', shell=True)
            time.sleep(0.1) 
elif choice == '2':
    for nxy in range(initial_nxy, final_nxy + 1):
            Data(-111, mcs, temperature, nxy, nxy, spins_orientation).write_params()
            subprocess.run(f'nohup ./xgenerate > ./generate_output/n{nxy}generate.output &', shell=True)
            time.sleep(0.1) 
