#oriented at https://www.usm.uni-muenchen.de/~ildar/hydro_course/ws2223/T06/setup_slab.pro by Klaus Dolag
import g3read as g3
import numpy as np
from shutil import copyfile

gamma=5./3.
rho=0.125
P = 1.0

output_file = 'wave.ic'
input_file = 'glass_10x10x10'

glass = g3.GadgetFile(input_file)
box_len = glass.header.BoxSize
nbase  = glass.header.npart[0]
grid = g3.read_new(input_file,'POS ', 0)

#we use 200 cubes of the base grid
number_of_boxes = 20
npart = number_of_boxes*nbase

pos = np.zeros((npart, 3), np.float32)
vel = np.zeros((npart, 3), np.float32)
m   = np.zeros(npart, np.float32)
U   = np.zeros(npart, np.float32)
#set up the low density-part
for i in range(0,number_of_boxes):
    l = i*nbase
    pos[l:l+nbase, 0] = grid[:,0] + i
    pos[l:l+nbase, 1] = grid[:,1] 
    pos[l:l+nbase, 2] = grid[:,2] 

id = np.array([i for i in range(npart)])+1

U[:] = P/(gamma-1)/rho
    
m[:] =  (rho*box_len**3)/npart

du = 0.01*P/(gamma-1)/rho #??
dv = 0.1
l=2.0
for j in range(npart):
    if pos[j,0] < l:
        vel[j,0] = dv*np.sin(np.pi*pos[j,0]/l)
    elif pos[j,0] > number_of_boxes-l:
        vel[j,0] = -1*dv*np.sin(np.pi*(number_of_boxes - pos[j,0])/l)

#dump coordiantes to text file for debugging
#coord = open('coord.txt', "w")
#for i in range(npart):
#    coord.write(str(pos[i,0]) + '\t' + str(pos[i,1]) + '\t' + str(pos[i,2]) )
#coord.close()

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
f.write_block("U   ", 0, U,   filename=output_file)