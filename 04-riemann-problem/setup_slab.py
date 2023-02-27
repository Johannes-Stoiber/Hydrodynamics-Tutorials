#oriented at https://www.usm.uni-muenchen.de/~ildar/hydro_course/ws2223/T04/slabsetup.f90 by Tadziu Hoffmann 
import g3read as g3
import numpy as np
from shutil import copyfile

#we use 10 cubes of the base grid for the low-density region and 8*reduced cubes of the base grid for the high density region
nbase = 1000
npart = 10*nbase+80*nbase

rho=1.5e-8
P=1.0, 
gamma=5.0/3.0
boxlen=1000.0

output_file = "slab.ic"

#read in the base grid 
glass = open('glass.txt')
rows = glass.readlines()
grid = np.zeros((nbase,3))
for i in range(nbase):
    zi = rows[i].split()
    grid[i,0] = zi[0]
    grid[i,1] = zi[1]
    grid[i,2] = zi[2]

pos = np.zeros((npart, 3), np.float32)
vel = np.zeros((npart, 3), np.float32)
m   = np.ones(npart, np.float32)
U   = np.zeros(npart, np.float32)
#set up the low density-part
for i in range(0,10):
    l = i*nbase
    pos[l:l+nbase, 0] = grid[:,0] + i
    pos[l:l+nbase, 1] = grid[:,1] 
    pos[l:l+nbase, 2] = grid[:,2] 
    
    vel[l:l+nbase, 0] = 0
    vel[l:l+nbase, 1] = 0
    vel[l:l+nbase, 2] = 0
    
    U[l:l+nbase] = 1.0
    
    m[l:l+nbase] =  1.0
    

#set up the high-density part
for i in range(0,20):
    for j in range(0,2):
        for k in range(0,2):
            l = (10+4*(i)+2*(j)+k)*nbase
            pos[l:l+nbase, 0] = 0.5*grid[:,0] + 0.5*(20+i)
            pos[l:l+nbase, 1] = 0.5*grid[:,1] + 0.5*(j)
            pos[l:l+nbase, 2] = 0.5*grid[:,2] + 0.5*(k)
            
            vel[l:l+nbase, 0] = 0
            vel[l:l+nbase, 1] = 0
            vel[l:l+nbase, 2] = 0
            
            U[l:l+nbase] = 1.0
    
            m[l:l+nbase] =  1.0

id = np.array([i for i in range(npart)])+1

#dump coordiantes to text file for debugging
#coord = open('coord.txt', "w")
#for i in range(npart):
#    coord.write(str(pos[i,0]) + '\t' + str(pos[i,1]) + '\t' + str(pos[i,2]) )
#coord.close()

copyfile("ic_header", output_file)

# read reference file    
f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
f.header.npart = [npart, 0, 0, 0, 0, 0]
f.header.npartTotal = [npart, 0, 0, 0, 0, 0]
#f.header.time = 0.0
f.header.boxsize = 1.0

f.write_header(f.header,filename=output_file)
    
f.add_file_block('POS ', (npart)*4*3  , partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('VEL ', (npart)*4*3  , partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('ID  ', (npart)*4 , partlen=4, dtype=np.int32) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 int)
f.add_file_block('MASS', (npart)*4 , partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)
f.add_file_block('U   ', (npart)*4 , partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)

# print(N**3, 'len pos',len(pos), 'len vel',len(vel), 'len IDs',len(id))
f.write_block("POS ", 0, pos, filename=output_file)
f.write_block("VEL ", 0, vel, filename=output_file)
f.write_block("ID  ", 0, id,  filename=output_file)
f.write_block("MASS", 0, m,   filename=output_file)
f.write_block("U   ", 0, U,   filename=output_file)
