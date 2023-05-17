import subprocess
import time 
import numpy as np

#Constants
g = 3.8
k = 0
iT = 1
dT = 0.1
fT = 100
order = 2
mag_dir = 'easy'
iH = 0
dH = 0.1
mcs = 5000
neq = 5000



seeds = [-48703, -23670, -13321, -67541, -91793]
print('Seeds: ', seeds)

print('Select the program:')
print('sLoop: (s)')
print('tLoop: (t)')
print('hLoop: (h)')
program = input('Else for exit:   ')


print('Please select a compound:')
compound = input('[1] <-- CrI3 ||| CrI2 --> [2]')
if compound == '1':
    compound = 'CrI3'
    J = 9.73
elif compound == '2':
    compound = 'CrI2'
    J = -4.73
else: 
    print('Invalid input, exiting...')
    exit()

print('Energies:')
print('J = ', J)
use_k = input('Use k?: (y/else)')
if use_k == 'y':
    k = 0.67
print('k = ', k)


initial_n = int(input('Enter the initial nx and ny (nx = ny): '))
final_n = int(input('Enter the final nx and ny (nx = ny): '))
step_n = int(input('Enter the step: '))

mcs = int(input('Enter the number of mcs: '))


def output_file(program, iseed, n, ):
    if program == 's':
        return f'{compound}n{n}t{iT}o{order}md{mag_dir}k{k}H{iH}s{iseed + 1}'
    elif program == 't':
        return f'{compound}n{n}dT{dT}o{order}md{mag_dir}k{k}H{iH}s{iseed + 1}'
    elif program == 'h':
        return f'{compound}n{n}dH{dH}o{order}md{mag_dir}k{k}t{iT}s{iseed + 1}'


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


def executable(program):
    if program == 's':
        return './xsLoop'
    elif program == 't':
        return './xtLoop'
    elif program == 'h':
        return './xhLoop'
    

number_of_simulations = len(seeds) * len(range(initial_n, final_n + 1, step_n))
print('You are going to print the following values for N:')
print(list(range(initial_n, final_n + 1, step_n)))
print(' also you are going to print the following values for seeds:')
print(seeds)
print('That would call a total of ', number_of_simulations, ' simulations.')
print('Are you sure you want to continue?')
print('Enter 1 to continue or anything else to exit: ')
confirmation = input()

if confirmation == '1':
    for n in range(initial_n, final_n + 1, step_n):
        for iseed in range(len(seeds)):
            print('Running simulation for: ', seeds[iseed], ' and N = ', n)
            write_input(program, iseed, n)
            exe = executable(program)
            log = f'{program}{output_file(program, iseed, n)}.log'
            subprocess.run(f'nohup {exe} > ../logs/{log} &', shell=True)
            time.sleep(0.1)
else: 
    print('Exiting...')
    exit()
        