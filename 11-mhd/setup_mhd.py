#oriented at https://www.usm.uni-muenchen.de/~ildar/hydro_course/ws2223/T04/slabsetup.f90 by Tadziu Hoffmann 
import g3read as g3
import numpy as np
from shutil import copyfile

input_file = 'glass_10x10x10'
glass = g3.GadgetFile(input_file)
h = glass.header
grid = g3.read_new(input_file,'POS ', 0)
#hglass = g3.read_new(input_file, 'HSML', 0)
boxsize = h.BoxSize
grid = grid/boxsize
nbase  = int(h.npart[0])

#we use 10 cubes of the base grid for the low-density region and 8*reduced cubes of the base grid for the high density region
nbase = 1000
npart = 50*nbase+400*nbase

output_file = "mhd.ic"


Psetupl=1.0
rhosetupl=1.0
rhosetupr=0.125
Psetupr=0.1

gamma=5.0/3.0
boxlen=100.0

pos = np.zeros((npart, 3), np.float32)
vel = np.zeros((npart, 3), np.float32)
m   = np.zeros(npart, np.float32)
U   = np.zeros(npart, np.float32)
b = np.zeros((npart, 3), np.float32)
hsml   = np.zeros(npart, np.float32)
rho   = np.zeros(npart, np.float32)
#set up the low density-part
for i in range(0,50):
    l = i*nbase
    pos[l:l+nbase, 0] = grid[:,0] + i
    pos[l:l+nbase, 1] = grid[:,1] 
    pos[l:l+nbase, 2] = grid[:,2] 
    
    vel[l:l+nbase, 0] = 0
    vel[l:l+nbase, 1] = 0
    vel[l:l+nbase, 2] = 0
    
    b[l:l+nbase, 0] = 0.75
    b[l:l+nbase, 1] = -1.0
    b[l:l+nbase, 2] = 0
    
    U[l:l+nbase] = Psetupr/(gamma-1)/rhosetupr
    
    #hsml[l:l+nbase]=hglass[i]
    
    #rho[l:l+nbase]=rhosetupr
    
    m[l:l+nbase] = rhosetupl/(2**3*nbase)



#set up the high-density part
for i in range(0,100):
    for j in range(0,2):
        for k in range(0,2):
            l = (50+4*(i)+2*(j)+k)*nbase
            pos[l:l+nbase, 0] = 0.5*grid[:,0] + 0.5*(100+i)
            pos[l:l+nbase, 1] = 0.5*grid[:,1] + 0.5*(j)
            pos[l:l+nbase, 2] = 0.5*grid[:,2] + 0.5*(k)
            
            vel[l:l+nbase, 0] = 0
            vel[l:l+nbase, 1] = 0
            vel[l:l+nbase, 2] = 0
            
            b[l:l+nbase, 0] = 0.75
            b[l:l+nbase, 1] = 1.0
            b[l:l+nbase, 2] = 0
            
            U[l:l+nbase] = Psetupl/(gamma-1)/rhosetupl
    
            #hsml[l:l+nbase]=hglass[i]
        
            #rho[l:l+nbase]=rhosetupr
            
            m[l:l+nbase] = rhosetupl/(2**3*nbase)
        
pos[:,0] -=50
pos[:,0] = -pos[:,0]
pos[:,0] +=50


id = np.array([i for i in range(npart)])+1


copyfile(r"C:\Users\johan\Python\Hydrodynamics\Tutorial01\ic_header", output_file)

# read reference file    
f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
f.header.npart = [npart, 0, 0, 0, 0, 0]
f.header.npartTotal = [npart, 0, 0, 0, 0, 0]
f.header.boxsize = 1

f.write_header(f.header,filename=output_file)
    
f.add_file_block('POS ', npart*4*3  , partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('VEL ', npart*4*3, partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('BFLD', npart*4*3, partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('ID  ', npart*4, partlen=4, dtype=np.int32) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 int)
f.add_file_block('MASS', npart*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)
f.add_file_block('U   ', npart*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)

# print(N**3, 'len pos',len(pos), 'len vel',len(vel), 'len IDs',len(id))
f.write_block("POS ", 0, pos, filename=output_file)
f.write_block("VEL ", 0, vel, filename=output_file)
f.write_block("BFLD", 0, b, filename=output_file)
f.write_block("ID  ", 0, id,  filename=output_file)
f.write_block("MASS", 0, m,   filename=output_file)
f.write_block("U   ", 0, U,   filename=output_file)