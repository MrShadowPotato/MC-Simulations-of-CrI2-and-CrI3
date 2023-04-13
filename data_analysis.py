seeds = [-111, -222, -333]
print('Seeds: ', seeds)

energies = []
magnetizations = []
print('The steps between ns are 5')
initial_n = int(input('Enter the initial n: '))
final_n = int(input('Enter the final n: '))
for n in range(initial_n, final_n + 1, 5):
    energy_sum = 0
    magnetization_sum = 0
    for seed in seeds:
        stabilizer_format = f'data/stabilizer/nxy{n}t5m3000s{seed}o1.txt'
        with open(stabilizer_format, 'r') as data:
            last_line = data.readlines()[-2]
            energy = float(last_line.split()[1])
            magnetization = float(last_line.split()[2])
        energy_sum += energy
        magnetization_sum += magnetization
    energy_avg = energy_sum/len(seeds)
    magnetization_avg = magnetization_sum/len(seeds)
    energies.append(energy_avg)
    magnetizations.append(magnetization_avg)
    open('em-size.txt', 'a').write(f'{n} {energy_avg} {magnetization_avg}')

        
