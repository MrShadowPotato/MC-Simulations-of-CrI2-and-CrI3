class Data:
    def __init__(self, seed, mcs, temperature, nx, ny, spins_orientation, output_file):
        self.seed = seed
        self.mcs = mcs
        self.temperature = temperature
        self.nx = nx
        self.ny = ny
        self.spins_orientation = spins_orientation
        self.output_file = output_file
    

with open ('params.txt', 'w') as f:
    pass

