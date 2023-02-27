import g3read as g3
import numpy as np
from shutil import copyfile

########read file##########
input_file = 'glass_10x10x10'
glass = g3.GadgetFile(input_file)
box_size = glass.header.BoxSize
nbase = glass.header.npart[0]
grid = g3.read_new(input_file,'POS ', 0) #works now
grid = grid/box_size

#####create output file#######

# free parameters
output_file = 'blast_100x100x100.ic'
nboxes = 10
nblast = 20

# creating...
npart = nboxes**3 * nbase

pos = np.zeros((npart, 3))

i = 0
for ix in range(nboxes):
    for iy in range(nboxes):
        for iz in range(nboxes):
            for j in range(nbase):
                pos[i,0] = (grid[j,0] + ix) / nboxes
                pos[i,1] = (grid[j,1] + iy) / nboxes
                pos[i,2] = (grid[j,2] + iz) / nboxes
                i += 1
                
pos = pos/(box_size)
#pos = pos - 0.5
#pos = pos*2
                
vel = np.zeros((npart, 3))
u = np.zeros(npart)
m = np.zeros(npart)
id = np.array([i for i in range(npart)])+1

gamma = 5./3.
rho = 0.125
Pback = 1e-20
Pblast = 1

u[:] = Pback/(gamma-1)/rho
m[:] = (rho*box_size**3)/npart

radii = np.sqrt( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 )
sradii = np.sort( radii )
maxr = sradii[nblast]
cond = np.where( (pos[:,0] - 0.5)**2 + (pos[:,1]-0.5)**2 + (pos[:,2]-0.5)**2 < maxr**2 )
u[cond] = Pblast/(gamma-1)/rho

copyfile('ic_header', output_file)

# read reference file    
f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
f.header.npart = [npart, 0, 0, 0, 0, 0]
f.header.npartTotal = [npart, 0, 0, 0, 0, 0]
f.header.time = 0.0
f.header.boxsize = 1.0

f.write_header(f.header,filename=output_file)
    
f.add_file_block('POS ', npart*4*3  , partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('VEL ', npart*4*3, partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('ID  ', npart*4, partlen=4, dtype=np.int32) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 int)
f.add_file_block('MASS', npart*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)
f.add_file_block('U   ', npart*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)

# print(N**3, 'len pos',len(pos), 'len vel',len(vel), 'len IDs',len(id))
f.write_block("POS ", 0, pos, filename=output_file)
f.write_block("VEL ", 0, vel, filename=output_file)
f.write_block("ID  ", 0, id,  filename=output_file)
f.write_block("MASS", 0, m,   filename=output_file)
f.write_block("U   ", 0, u,   filename=output_file)      