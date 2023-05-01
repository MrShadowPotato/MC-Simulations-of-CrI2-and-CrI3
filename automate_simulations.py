import subprocess
import time

class Data:
    def __init__(self, seed = -111, mcs = 5000, temperature = 5, nx = 0, ny = 0, spins_orientation = 1, iT = 1, fT=100, dT=1, idx=1000, compound='CrI3',md='easy'):
        self.compound = compound
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
        self.magnetization_direction = md
        self.path_simulation_data = self.path_simulation_data()
        self.path_temperature_data = self.path_temperature_data()

    def path_simulation_data(self):
        if self.spins_orientation == '1':
            return f'{self.compound}n{self.nx}t{self.temperature}m{self.mcs}o{self.spins_orientation}s{self.seed}.txt'
        else:
            return f'{self.compound}n{self.nx}t{self.temperature}m{self.mcs}o{self.spins_orientation}{self.magnetization_direction}s{self.seed}.txt'
    
    def path_temperature_data(self):
        return f'{self.compound}n{self.nx}iT{self.iT}fT{self.fT}dT{self.dT}m{self.mcs}i{self.idx}o{self.spins_orientation}s{self.seed}.txt'
    
    def write_params(self):
        with open ('params.txt', 'w') as f:
            f.write(f'seed = {self.seed}\n')
            f.write(f'mcs = {self.mcs}\n')
            f.write(f'temperature = {self.temperature}\n')
            f.write(f'nx = {self.nx}\n')
            f.write(f'ny = {self.ny}\n')
            f.write(f'spins_orientation = {self.spins_orientation}\n')
            f.write(f'path_simulation_data = {self.path_simulation_data}\n')
            f.write(f'Cr_xyz = n{self.nx}Cr_{self.compound}.xyz\n')
            f.write(f'Cr_neighbors = n{self.nx}Cr_{self.compound}neighbors.txt\n')
            f.write(f'initial_temperatre = {self.iT}\n')
            f.write(f'final_temperature = {self.fT}\n')
            f.write(f'temperature_step = {self.dT}\n')
            f.write(f'index = {self.idx}\n')
            f.write(f'path_temperature_data = {self.path_temperature_data}\n')
            f.write(f'compound = {self.compound}\n')
            f.write(f'magnetization_direction = {self.magnetization_direction}')


print('Please select a compound:')
compound = input('[1] <-- CrI3 ||| CrI2 --> [2]')
if compound == '1':
    compound = 'CrI3'
elif compound == '2':
    compound = 'CrI2'
else: 
    print('Invalid input, exiting...')
    exit()

initial_nxy = int(input('Enter the initial nx and ny (nx = ny): '))
final_nxy = int(input('Enter the final nx and ny (nx = ny): '))
step_n = int(input('Enter the step: '))


#seeds = [-111, -222, -333, -444, -555]
seeds = [-111]
print('Seeds: ', seeds)
print('Enter 1 to run xgenerate')

print('Enter 2 to run xsimulation')
print('Enter 3 to run xtemp')

choice = input('Choice:   ')


#Program selection
if choice == '1':   #xgenerate
    pass
elif choice == '2' or choice == '3':
    
    if choice == '2':   #xsimulation
        print('Now choose what hamiltonian will you use:')
        print('Only exchange\n--------> [1]\n')
        print('Exchange + anisotrpy\n--------> [2]\n')
        hamiltonian = input('Enter the hamiltonian:   ')
        if hamiltonian == '1':
            program = 'stabilizer'
        elif hamiltonian == '2':
            program = 'stabilizer_anisotropy'
        else:
            print('Invalid input, exiting...')
            exit()
        temperature = int(input('Enter the temperature:   '))
        print('Now choose the orientation for initial magnetization:  ')
        spins_orientation = int(input('1 for all random and 2 for all the same:  '))
    elif choice == '3':   #xtemp
        iT = int(input('Enter the initial temperature:   '))
        fT = int(input('Enter the final temperature:   '))
        dT = int(input('Enter the temperature step:   '))
        idx = int(input('Enter the index averager:   '))
        spins_orientation = 2

    if spins_orientation == 1:
        print('You have chosen all random spins')
    elif spins_orientation == 2:
        print('\n'*2,'*'*20,'\n', 'Now choose the direction for initial magnetization: ')
        print('For Easy axis direction press\n--------> [1]\n')
        print('For Reversed Easy axis direction press\n--------> [2]\n')
        print('For x (p1) axis direction press\n--------> [3]\n')
        print('For y (p2) axis direction press\n--------> [4]\n')
        print('For xy plane (p1 + p2) direction press\n--------> [5]\n')
        print('For 45deg between xy plance and easy axis direction press\n--------> [6]\n')

        magnetization_choice = input()

        if magnetization_choice == '1':
            magnetization_direction = 'easy'
        elif magnetization_choice == '2':
            magnetization_direction = 'reversed_easy'
        elif magnetization_choice == '3':
            magnetization_direction = 'p1'
        elif magnetization_choice == '4':
            magnetization_direction = 'p2'
        elif magnetization_choice == '5':
            magnetization_direction = 'p1_p2plane'
        elif magnetization_choice == '6':
            magnetization_direction = '45deg'
        else:
            print('Invalid input, exiting...')
            exit()    

    mcs = int(input('Enter the number of Monte Carlo steps:  '))      
else:
    print('Invalid input, exiting...')
    exit()


number_of_simulations = len(seeds)*len(list(range(initial_nxy, final_nxy + 1, step_n)))

print('You are going to print the following values for N:')
print(list(range(initial_nxy, final_nxy + 1, step_n)))
print(' and you are going to print the following values for seeds:')
print(seeds)
print('That would call a total of ', number_of_simulations, ' simulations.')
print('Simulations : ', number_of_simulations)
print('Are you sure you want to continue?')
print('Enter 1 to continue or anything else to exit: ')
confirmation = input()


if confirmation == '1':
    for nxy in range(initial_nxy, final_nxy + 1, step_n):
        if choice == '1':
            Data(nx=nxy, ny=nxy, compound=compound).write_params()
            subprocess.run(f'nohup ./xgenerate > ./logs/structures/n{nxy}generate{compound}.output &', shell=True)
            time.sleep(0.1) 
        elif choice == '2' or choice == '3':  
            print(f'Running for nxy = {nxy}...')
            for seed in seeds:
                print(f'Running for seed = {seed}...')
                if choice == '2':   #Stabilizer
                    if spins_orientation == 1:
                        Data(seed=seed, mcs=mcs, temperature=temperature, nx=nxy, ny=nxy, spins_orientation=spins_orientation, compound=compound).write_params()
                    elif spins_orientation == 2:
                        Data(seed=seed, mcs=mcs, temperature=temperature, nx=nxy, ny=nxy, spins_orientation=spins_orientation, compound=compound, md=magnetization_direction).write_params()
                    subprocess.run(f'nohup ./x{program} > ./logs/{program}/n{nxy}t{temperature}s{seed}.output &', shell=True)
                elif choice == '3':  #Temperature
                    Data(spins_orientation=spins_orientation, iT=iT, fT=fT, dT=dT, idx=idx, mcs=mcs, nx=nxy, ny=nxy, seed=seed, compound=compound, md=magnetization_direction).write_params()
                    subprocess.run(f'nohup ./xtemp > ./logs/temperature/n{nxy}t{iT}-{fT}d{dT}s{seed}.output &', shell=True)
                time.sleep(0.1) 
                        
else: 
    print('Exiting...')
    exit()
