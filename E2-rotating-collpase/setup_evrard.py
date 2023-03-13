import g3read as g3
import numpy as np
from shutil import copyfile

M_star=1.5815
R_star=1.0
G=1.0

########read file##########
input_file = 'glass_100x100x100'
glass = g3.GadgetFile(input_file)
box_size = glass.header.BoxSize
Nin = glass.header.npart[0]
pos = glass.read_new('POS ', 0)
#######end read file#########

#####create output file#######
output_file = 'evrard_shear.ic'
N=20000


               
vel = np.zeros((N,3))
u = np.zeros(N)
m = np.zeros(N)
id = np.array([i for i in range(N)])+1

####ensure everything is within [-1,1]                
pos = pos/box_size
pos = pos - 0.5
pos = pos*2
 
####compute the radii###
radii = np.sqrt( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 )
                      
####scale coordinates with radius to create r profile
for i in range(0,3):
    pos[:,i] = pos[:,i]*radii
    #pos[:,i] = np.exp(pos[:,i])
                      
###we want to have N particles in a sphere
rr = np.sqrt( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 )
srr = np.sort(rr)
max_r = srr[N]
cond = np.where( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 < max_r**2 )
pos = pos[cond]

rr_again = np.sqrt( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 )
##rescale to the above units
pos = pos/np.max(rr_again)*R_star
m[:]=M_star/N
u[:]=G*M_star/R_star*0.05


# Proper rotating attempt now
velocity_amplitude = 1.0
yx_radii = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
thetas = np.arctan2(pos[:,1], pos[:,0])
cyl2cart = np.array([-1*np.sin(thetas) , np.cos(thetas), np.zeros(pos[:,2].shape)])
rotating_velocities = velocity_amplitude * (yx_radii * cyl2cart)

vel = rotating_velocities.T

# Create shear
sign_changer = np.array([pos[:,2]/np.abs(pos[:,2])]).T # Column of ones where z>0 and negative ones where z<0)
vel = sign_changer * vel # Switches velocities to negative current value where z<0

copyfile('ic_header', output_file)

# read reference file    
f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
f.header.npart = [N, 0, 0, 0, 0, 0]
f.header.npartTotal = [N, 0, 0, 0, 0, 0]
f.header.boxsize = 1.0

f.write_header(f.header,filename=output_file)
    
f.add_file_block('POS ', N*4*3  , partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('VEL ', N*4*3, partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
f.add_file_block('ID  ', N*4, partlen=4, dtype=np.int32) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 int)
f.add_file_block('MASS', N*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)
f.add_file_block('U   ', N*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)

# print(N**3, 'len pos',len(pos), 'len vel',len(vel), 'len IDs',len(id))
f.write_block("POS ", 0, pos, filename=output_file)
f.write_block("VEL ", 0, vel, filename=output_file)
f.write_block("ID  ", 0, id,  filename=output_file)
f.write_block("MASS", 0, m,   filename=output_file)
f.write_block("U   ", 0, u,   filename=output_file)  
