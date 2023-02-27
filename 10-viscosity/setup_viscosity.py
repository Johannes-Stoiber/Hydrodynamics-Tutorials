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
nboxes_x = 10
nboxes_y = 30
nboxes_z = 1

output_file = 'viscosity.ic'

npart = nboxes_x*nboxes_y*nboxes_z*nbase

pos = np.zeros((npart, 3), np.float32)

# setup initial particles
i = 0
for ix in range(nboxes_x):
    for iy in range(nboxes_y):
        for iz in range(nboxes_z):
            for j in range(nbase):
                pos[i,0] = (grid[j,0] + ix) / nboxes_z
                pos[i,1] = (grid[j,1] + iy) / nboxes_z
                pos[i,2] = (grid[j,2] + iz) / nboxes_z
                i += 1


vel = np.zeros((npart, 3), np.float32)
u = np.zeros(npart)
m = np.zeros(npart)           
id = np.array([i for i in range(npart)])+1

gamma = 5./3.
rho = 1.0
P = 1.0

u_ini = P/(gamma-1)/rho
u[:] = u_ini
mpart = (rho * nboxes_x**3 * boxsize**3)/npart
m[:] = mpart

# setup bottom wall
b = np.where(pos[:,1] < 1)
id[b] = 0
m[b] = mpart*4
u[b] = 1.1*u_ini

#decrese mass of top particles
top = np.where(pos[:,1] > 10)
m[top] = mpart/5
u[top] = 5*u_ini

#create tube
t1 = np.where((pos[:,1] > 12) & (pos[:,0] > 3)  & (pos[:,0] < 4))
id[t1] = 0
m[t1] = mpart*4
u[t1] = 1.1*u_ini
t2 = np.where((pos[:,1] > 12) & (pos[:,0] > 6)  & (pos[:,0] < 7))
id[t2] = 0
m[t2] = mpart*4
u[t2] = 1.1*u_ini

#create moving wall
dv = -0.1
w = np.where((pos[:,1] > 29) & (pos[:,0] > 4) & (pos[:,0] < 6))
id[w] = 0
vel[:,1][w] = dv
m[w] = 8*mpart
u[w] = 2*u_ini

#fluid within tube
f = np.where( (pos[:,1] > 12) & (pos[:,1] < 29)  & (pos[:,0] > 4) & (pos[:,0] < 6) )
m[f] = 5*mpart
vel[:,1][f] = dv
u[f] = u_ini/5

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