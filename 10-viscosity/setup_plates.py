#oriented at https://www.usm.uni-muenchen.de/~ildar/hydro_course/ws2223/T08/setup_nozzle.pro by Klaus Dolag 
import g3read as g3
import numpy as np
from shutil import copyfile

input_file = 'glass_10x10x10'
glass = g3.GadgetFile(input_file)
h = glass.header
grid = g3.read_new(input_file,'POS ', 0)
boxsize = h.BoxSize
grid = grid/boxsize
nbase  = int(h.npart[0])

#free parameters
nboxes_xy = 10
nboxes_z = 1

output_file = 'plates.ic'

npart = nboxes_xy**2*nboxes_z*nbase

pos = np.zeros((npart, 3), np.float32)

# setup initial particles
i = 0
for ix in range(nboxes_xy):
    for iy in range(nboxes_xy):
        for iz in range(nboxes_z):
            for j in range(nbase):
                pos[i,0] = (grid[j,0] + ix) / nboxes_xy 
                pos[i,1] = (grid[j,1] + iy) / nboxes_xy
                pos[i,2] = (grid[j,2] + iz) / nboxes_xy
                i += 1


vel = np.zeros((npart, 3), np.float32)
u = np.zeros(npart)
m = np.zeros(npart)           
id = np.array([i for i in range(npart)])+1

gamma = 5./3.
rho = 1.0
P = 1.0
Pw = 1.05

u[:] = P/(gamma-1)/rho
mpart = (rho * nboxes_xy**3 * boxsize**3)/npart
m[:] = mpart

# setup plates
dv = 0.5
L = 0.1
p1 = np.where(pos[:,0] < L)
p2 = np.where(pos[:,0] > (1 - L))
id[p1] = 0
id[p2] = 0
vel[:,1][p1] = dv
vel[:,1][p2] = -dv
u[p1] = Pw/(gamma-1)/rho
u[p2] = Pw/(gamma-1)/rho


copyfile(r"ic_header", output_file)

# read reference file    
f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
f.header.npart = [npart, 0, 0, 0, 0, 0]
f.header.npartTotal = [npart, 0, 0, 0, 0, 0]
f.header.boxsize = 1

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