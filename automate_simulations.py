import subprocess
import time

class Data:
    def __init__(self, seed = -111, mcs = 5000, temperature = 5, nx = 0, ny = 0, spins_orientation = 1, iT = 1, fT=100, dT=1, idx=1000):
        self.seed = seed
        self.mcs = mcs
        self.temperature = temperature
        self.nx = nx
        self.ny = ny
        self.spins_orientation = spins_orientation
        self.iT = iT
        self.fT = fT
        self.dT = dT
        self.idx = idx
        self.output_file = self.output_file()
        self.temp_iterator_file = self.temp_iterator_file()

    def output_file(self):
        return f'nxy{self.nx}t{self.temperature}m{self.mcs}s{self.seed}o{self.spins_orientation}.txt'
    
    def temp_iterator_file(self):
        return f'n{self.nx}t{self.iT}-{self.fT}d{self.dT}m{self.mcs}i{self.idx}o{self.spins_orientation}s{self.seed}.txt'
    
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
            f.write(f'initial_temperatre = {self.iT}\n')
            f.write(f'final_temperature = {self.fT}\n')
            f.write(f'temperature_step = {self.dT}\n')
            f.write(f'index = {self.idx}\n')
            f.write(f'temp_iterator_file = {self.temp_iterator_file}\n')
    



#seeds = [-111, -222, -333, -444, -555]
seeds = [-111]
print('Seeds: ', seeds)
print('Enter 1 to run xsimulation')
print('Enter 2 to run xgenerate')
print('Enter 3 to run xtemp')
choice =input('Choice: ')


initial_nxy = int(input('Enter the initial nx and ny (nx = ny): '))
final_nxy = int(input('Enter the final nx and ny (nx = ny): '))
step_n = int(input('Enter the step: '))


#Inputs
if choice == '1':
    mcs = int(input('Enter the number of Monte Carlo steps: '))        
    temperature = int(input('Enter the temperature: '))
    print('Enter the spins orientation: ')
    spins_orientation = int(input('1 for all random and 2 for all the same:'))

elif choice == '3':
    mcs = int(input('Enter the number of Monte Carlo steps: '))
    spins_orientation = 2
    iT = int(input('Enter the initial temperature: '))
    fT = int(input('Enter the final temperature: '))
    dT = int(input('Enter the temperature step: '))
    idx = int(input('Enter the index averager: '))

number_of_simulations = len(seeds)*len(list(range(initial_nxy, final_nxy + 1, step_n)))
print('You are going to print the following values for N:')
print(list(range(initial_nxy, final_nxy + 1, step_n)))
print(' and you are going to print the following values for seeds:')
print(seeds)
print('That would call a total of ', number_of_simulations, ' simulations.')
print('Simulations : ', number_of_simulations)
print('Are you sure you want to continue?')
print('Enter 1 to continue or anything else to exit: ')

confirmation = int(input())

if confirmation == 1:
    if choice == '1' or choice == '3':
            for nxy in range(initial_nxy, final_nxy + 1, step_n):
                print(f'Running for nxy = {nxy}...')
                for seed in seeds:
                    print(f'Running for seed = {seed}...')
                    if choice == '1':
                        Data(seed, mcs, temperature, nxy, nxy, spins_orientation).write_params()
                        subprocess.run(f'nohup ./xsimulation > ./stabilizer_output/n{nxy}t{temperature}s{seed}.output &', shell=True)
                    elif choice == '3':
                        Data(spins_orientation=spins_orientation, iT=iT, fT=fT, dT=dT, idx=idx, mcs=mcs, nx=nxy, ny=nxy, seed=seed).write_params()
                        subprocess.run(f'nohup ./xtemp > ./temp_output/n{nxy}t{iT}-{fT}d{dT}s{seed}.output &', shell=True)
                    time.sleep(0.1) 
                        
    elif choice == '2':
        for nxy in range(initial_nxy, final_nxy + 1):
                Data(-111, mcs, temperature, nxy, nxy, spins_orientation).write_params()
                subprocess.run(f'nohup ./xgenerate > ./generate_output/n{nxy}generate.output &', shell=True)
                time.sleep(0.1) 


else: 
    print('Exiting...')
    exit()
