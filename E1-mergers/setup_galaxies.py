import g3read as g3
import numpy as np

npart = np.array([15238, 141380, 60952, 32610, 0, 0], dtype = np.int32)
npart2 = 2*npart

data = np.loadtxt('galaxy.txt', dtype = np.float32)
ids0 = data[:,0].astype(np.int32)
m0 = data[:,1]
u0 = data[:,2]
pos0 = data[:,3:6]
vel0 = data[:,6:]

nbase = np.sum(npart)

ntot = 2*nbase

pos = np.zeros((ntot, 3))
vel = np.zeros((ntot, 3))
m = np.zeros(ntot)
u = np.zeros(ntot)

#rot_counter = np.array([[1,0,0],[0,-1,0],[0,0,-1]])
#rot_same = np.array([[1,0,0],[0,1,0],[0,0,1]])
#rot = np.array([[1,0,0],[0,-1,0],[0,0,-1]])
rotless = np.array([[1,0,0],[0,-0.99,0.09],[0,-0.09,-0.99]])

i = 0
j = 0
for l in range(4):
    n = npart[l]
    n2 = npart2[l]
    #print('l =',l,'n =',n,'n2=',n2,'i=',i )
    pos[i:i+n] = pos0[j:j+n] + np.array([0,-50,0])
    pos[i+n:i+n2] = np.matmul(pos0[j:j+n],rotless) + np.array([0,50,0])
    vel[i:i+n] = vel0[j:j+n] + np.array([0,10,0])
    vel[i+n:i+n2] = np.matmul(vel0[j:j+n],rotless) + np.array([0,-10,0]) 
    m[i:i+n] = m0[j:j+n]
    m[i+n:i+n2] = m0[j:j+n]
    u[i:i+n] = u0[j:j+n]
    u[i+n:i+n2] = u0[j:j+n]  
    i += n2
    j += n
    
ids = np.array([i for i in range(ntot)])+1

output_file = 'galaxies_counter.ic'

with open(output_file, "w") as f:
    pass

redshift = 1.0
time = 0.0
Omega0 = 0.27
OmegaLambda = 1.0 - Omega0
HubbleParam = 0.704
num_files = 1
BoxSize = 1.0

mass_array = np.array([0, m[npart2[0]], m[npart2[0]+npart2[1]], m[npart2[0]+npart2[1]+npart2[2]],0,0], dtype=np.float32)

header = g3.GadgetHeader(npart2, mass_array, time, redshift, BoxSize, Omega0, OmegaLambda, HubbleParam, num_files=num_files)
f = g3.GadgetWriteFile(output_file, npart2, {}, header)
f.write_header(f.header)

f.add_file_block('ID  ', ntot * 4, partlen = 4, dtype=np.int32)
f.add_file_block('POS ', ntot * 4 * 3, partlen = 4 * 3)
f.add_file_block('VEL ', ntot * 4 * 3, partlen = 4 * 3)
f.add_file_block('MASS', npart2[0] * 4, partlen = 4)
f.add_file_block('U   ', npart2[0] * 4, partlen = 4)

f.write_block('ID  ', -1, ids)
f.write_block('POS ', -1, pos)
f.write_block('MASS', -1, m[:npart2[0]])
f.write_block('VEL ', -1, vel)
f.write_block('U   ', -1, u[:npart2[0]])