#oriented at https://www.usm.uni-muenchen.de/~ildar/hydro_course/ws2223/T04/slabsetup.f90 by Tadziu Hoffmann 
import g3read as g3
import numpy as np
from shutil import copyfile

gamma=5./3.
rho = 0.125
P = 1.0

input_file = 'glass_10x10x10'
glass = g3.GadgetFile(input_file)
box_len = glass.header.BoxSize
nbase  = glass.header.npart[0]
grid = g3.read_new(input_file,'POS ', 0)

#free parameters
output_file = 'nozzle.ic'
nboxes_x = 200
nboxes_yz = 2
npart = nboxes_x*nboxes_yz**2*nbase
dv = 2

pos = np.zeros((npart, 3), np.float32)
vel = np.zeros((npart, 3), np.float32)
m   = np.zeros(npart, np.float32)
u   = np.zeros(npart, np.float32)

# setup initial particles
i = 0
for ix in range(nboxes_x):
    for iy in range(nboxes_yz):
        for iz in range(nboxes_yz):
            for j in range(nbase):
                pos[i,0] = (grid[j,0] + ix) / nboxes_yz + 100
                pos[i,1] = (grid[j,1] + iy) / nboxes_yz
                pos[i,2] = (grid[j,2] + iz) / nboxes_yz
                i += 1
                
id = np.array([i for i in range(npart)])+1

# setup right wall
w1 = np.where(pos[:,0] > 195)
id[w1] = 0

# setup left wall
w2 = np.where( (pos[:,0] > 100) & (pos[:,0] < 105) )
id[w2] = 0
vel[w2] = -dv

# setup noozle
# setup noozle
l = 35
dn = 0.3
n1 = np.where( (pos[:,0] > 105) & (pos[:,0] < 140) & (pos[:,1] < dn*np.sin( np.pi*(pos[:,0]-105)/l ) ) ) 
n2 = np.where( (pos[:,0] > 105) & (pos[:,0] < 140) & (pos[:,1] > 1-dn*np.sin( np.pi*(pos[:,0]-105)/l ) ) ) 
id[n1] = 0
id[n2] = 0

pos[:,0] = pos[:,0]/5 #--> dv = 1

u[:] = P/(gamma-1)/rho
mpart = (rho * nboxes_yz^3 * h.BOXSIZE^3)/Np
m[:] = mpart

copyfile(r"C:\Users\johan\Python\Hydrodynamics\Tutorial01\ic_header", output_file)

# read reference file    
f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
f.header.npart = [npart, 0, 0, 0, 0, 0]
f.header.npartTotal = [npart, 0, 0, 0, 0, 0]
#f.header.time = 0.0
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