import g3read as g3
import numpy as np
from shutil import copyfile

########read file##########
#input_file = 'glass.txt'
input_file = 'glass_10x10x10'

glass = g3.GadgetFile("glass_10x10x10")
box_size = glass.header.BoxSize
Nin = glass.header.npart[0]
grid = g3.read_new(input_file,'POS ', 0) #somehow doesn't work

#read in the base grid 
#glass = open(input_file)
#rows = glass.readlines()
#grid = np.zeros((Nin,3))
#for i in range(Nin):
#    zi = rows[i].split()
#    grid[i,0] = zi[0]
#    grid[i,1] = zi[1]
#    grid[i,2] = zi[2]
    
grid = grid/box_size

#######end read file#########

#####create output file#######
output_file = 'glass_100x100x100'
Ndup = 10

Np = Ndup**3 * Nin

pos = np.zeros((Np, 3))

i = 0
for ix in range(Ndup):
    for iy in range(Ndup):
        for iz in range(Ndup):
            for j in range(Nin):
                pos[i,0] = (grid[j,0] + ix) / Ndup
                pos[i,1] = (grid[j,1] + iy) / Ndup
                pos[i,2] = (grid[j,2] + iz) / Ndup
                i += 1
                
vel = np.zeros((Np, 3))
u = np.zeros(Np)
m = np.zeros(Np)
id = np.array([i for i in range(Np)])+1

gamma = 5./3.
rho = 0.125
P = 1.0

u[:] = P/(gamma-1)/rho
m[:] = (rho*box_size**3)/Np


copyfile("ic_header", output_file)

# read reference file    
f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
f.header.npart = [Np, 0, 0, 0, 0, 0]
f.header.npartTotal = [Np, 0, 0, 0, 0, 0]
#f.header.time = 0.0
f.header.boxsize = 1.0

f.write_header(f.header,filename=output_file)
    
f.add_file_block('POS ', Np*4*3  , partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('VEL ', Np*4*3, partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('ID  ', Np*4, partlen=4, dtype=np.int32) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 int)
f.add_file_block('MASS', Np*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)
f.add_file_block('U   ', Np*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)

# print(N**3, 'len pos',len(pos), 'len vel',len(vel), 'len IDs',len(id))
f.write_block("POS ", 0, pos, filename=output_file)
f.write_block("VEL ", 0, vel, filename=output_file)
f.write_block("ID  ", 0, id,  filename=output_file)
f.write_block("MASS", 0, m,   filename=output_file)
f.write_block("U   ", 0, u,   filename=output_file)           
