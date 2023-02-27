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
nboxes_x = 40
nboxes_yz = 2
L = nboxes_x / nboxes_yz

print(f'Length of the region to fill is {L:.1}')

output_file = 'ac.ic'

npart = nboxes_x*nboxes_yz**2*nbase

pos = np.zeros((npart, 3), np.float32)

# setup initial particles
i = 0
for ix in range(nboxes_x):
    for iy in range(nboxes_yz):
        for iz in range(nboxes_yz):
            for j in range(nbase):
                pos[i,0] = (grid[j,0] + ix) / nboxes_yz 
                pos[i,1] = (grid[j,1] + iy) / nboxes_yz
                pos[i,2] = (grid[j,2] + iz) / nboxes_yz
                i += 1

vel = np.zeros((npart, 3), np.float32)
u = np.zeros(npart)
m = np.zeros(npart)           
id = np.array([i for i in range(npart)])+1

gamma = 5./3.
rho = 0.125
P = 1.0

u[:] = P/(gamma-1)/rho
mpart = (rho * nboxes_yz**3 * boxsize**3)/npart
m[:] = mpart

# setup first right wall
dL = 0.5
wr1 = np.where((pos[:,0] > (L-dL)) & (pos[:,1] > 0.66))
wr2 = np.where((pos[:,0] > (L-dL)) & (pos[:,1] < 0.33)) 
id[wr1] = 0
id[wr2] = 0

# setup second right wall
dL = 0.5
dv = 0.3
wr3 = np.where((pos[:,0] > (L-2*dL)) & (pos[:,0] < (L-dL)))
id[wr3] = 0
vel[wr3] = -dv

# setup nozzle 1
dL = 4
x0 = 7
dmin = 0.33
amp = (0.5-dmin)*2
dy = abs(pos[:,1]-0.5)
#dy_nozzle = 0.5 + amp * (np.sin((pos[:,0] - x0)/dL*2*np.pi + np.pi/2) - 1) 
n1 = np.where( (pos[:,0] > x0) & (pos[:,0] < x0+dL) & (dy > 0.5*amp) ) 
id[n1] = 0

# setup nozzle 2
dL = 4
x0 = 7
dmin = 0.33
dy_nozzle = dmin + 2*dmin/dL*(pos[:,0]-x0-dL)
dy = abs(pos[:,1]-0.5)
n1 = np.where( (pos[:,0] > x0+dL) & (pos[:,0] < x0+2*dL) & (dy > 0.5*dy_nozzle) ) 
id[n1] = 0

#move everything to the right by L
#pos[:,0] = pos[:,0] + L

copyfile(r"C:\Users\johan\Python\Hydrodynamics\Tutorial01\ic_header", output_file)

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