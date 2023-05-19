import subprocess
import time 
import numpy as np

#Constants
k = 0
iH = 0
iT = 1
g = 3.8
mcs = 5000
neq = 5000
order = 2
mag_dir = 'easy'
fT = 100
dT = 0.1
dH = 0.1
dS = 0.9
name_choice = 'no'



seeds = [-48703, -23670, -13321, -67541, -91793]
print('Seeds: ', seeds)

print('Select the program:')
print('sLoop:   (s)')
print('tLoop:   (t)')
print('hLoop:   (h)')
print('CrIxgen: (g)')
program = input('Else for exit:   ')


print('Please select a compound:')
compound = input('[1] <-- CrI3 ||| CrI2 --> [2]  ')
if compound == '1':
    compound = 'CrI3'
    J = 9.73
elif compound == '2':
    compound = 'CrI2'
    J = -4.73
else: 
    print('Invalid input, exiting...')
    exit()

initial_n = int(input('Enter the initial nx and ny (nx = ny): '))
final_n = int(input('Enter the final nx and ny (nx = ny): '))
step_n = int(input('Enter the step: '))

if program != 'g':
    print('Energies:')
    print('J = ', J)
    use_k = input('Use k?: (y/else)')
    if use_k == 'y':
        k = 0.67
    print('k = ', k)



    iT = float(input('Enter the iT (T): '))
    if program == 't':
        dT = float(input('Enter the dT: '))
        fT = float(input('Enter the fT: '))
    iH = float(input('Enter the iH (H): '))
    if program == 'h':
        dH = float(input('Enter the dH: '))



    mcs = int(input('Enter the number of mcs: '))
    if program == 'h' or program =='t':
        neq = int(input('Enter the number of equilibration steps: '))
        order = 2
    elif program == 's':
        print('Enter (1) for all random')
        print('Enter (2) for all the same')
        order = int(input('Enter the spin order: '))


    if order == 1:
        print('You have chosen all random spins')
    elif order == 2:
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
    else :
        print('Invalid input, exiting...')
        exit()



    print('For hLoop and tLoop dS equal to 0.2 and 0.9 respectively is recommended.')
    dS = input('Enter the dS (0,1] : ')



    print('Want to enter a custom output file name?')
    print('Enter (y) for yes or else for no')
    name_choice = input('Enter your choice: ')
    if name_choice == 'y':
        output_file_name = input('Enter the output file name: ')


def output_file(program, iseed, n):
    if name_choice == 'y':
        return f'{compound}n{n}{output_file_name}s{iseed + 1}' 
    elif program == 's':
        return f'{compound}n{n}t{iT}o{order}md{mag_dir}k{k}H{iH}s{iseed + 1}'
    elif program == 't':
        return f'{compound}n{n}dT{dT}o{order}md{mag_dir}k{k}H{iH}s{iseed + 1}'
    elif program == 'h':
        return f'{compound}n{n}dH{dH}o{order}md{mag_dir}k{k}t{iT}s{iseed + 1}'
    elif program == 'g':
        return f'n{n}{compound}'


def write_input(program, iseed, n):
    with open('input', 'w') as f:
        f.write('Name of the file where the output data will be stored:\n')
        output_file_name = output_file(program, iseed, n)
        f.write(f'{output_file_name}\n')
        f.write('**** System parameters ****\n')
        f.write(f'nx= {n}\n')
        f.write(f'ny= {n}\n')
        f.write(f'compound= {compound}\n')
        f.write(f'J= {J}\n')
        f.write(f'K= {k}\n')
        f.write(f'iH= {iH}\n')
        f.write(f'iT= {iT}\n')
        f.write(f'g= {g}\n')
        f.write('**** Montecarlo parameters ****\n')
        f.write(f'mcs= {mcs}\n')
        f.write(f'neq= {neq}\n')
        f.write(f'seed= {seeds[iseed]}\n')
        f.write('**** Initial spin config ****\n')
        f.write(f'order= {order}\n')
        f.write(f'mag_dir= {mag_dir}\n')
        f.write('******************************************\n')
        f.write('Data for tLoop \n')
        f.write(f'fT= {fT}\n')
        f.write(f'dT= {dT}\n')
        f.write('Data for hLoop \n')
        f.write(f'dH= {dH}\n')
        f.write('**** Spin multiplier ****\n')
        f.write(f'dS= {dS}')

def executable(program):
    if program == 's':
        return './xsLoop'
    elif program == 't':
        return './xtLoop'
    elif program == 'h':
        return './xhLoop'
    elif program == 'g':
        return './xgen'

number_of_simulations = len(range(initial_n, final_n + 1, step_n))
if program != 'g':
    number_of_simulations = len(seeds) * number_of_simulations
    print(' You are going to print the following values for seeds:')
    print(seeds)
print('You are going to print the following values for N:')
print(list(range(initial_n, final_n + 1, step_n)))
print('That would call a total of ', number_of_simulations, ' simulations.')
print('Are you sure you want to continue?')
print('Enter (y) to continue or anything else to exit: ')
confirmation = input()

if confirmation == 'y':
    for n in range(initial_n, final_n + 1, step_n):
        for iseed in range(len(seeds)):
            if program == 'g' and iseed != 0:
                continue                
            print('Running simulation for: ', seeds[iseed], ' and N = ', n)
            write_input(program, iseed, n)
            exe = executable(program)
            log = f'{program}{output_file(program, iseed, n)}.log'
            subprocess.run(f'nohup {exe} > ../logs/{log} &', shell=True)
            time.sleep(0.1)
else: 
    print('Exiting...')
    exit()
        