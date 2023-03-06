import g3read as g3
import numpy as np

npart = np.array([15238, 141380, 60952, 32610, 0, 0], dtype = np.int32)

data = np.loadtxt('galaxy.txt', dtype = np.float32)
id = data[:,0].astype(np.int32)
m = data[:,1]
u = data[:,2]
pos = data[:,3:6]
vel = data[:,6:]

ntot = np.sum(npart)

output_file = 'galaxy.ic'

with open(output_file, "w") as f:
    pass

redshift = 1.0
time = 0.0
Omega0 = 0.27
OmegaLambda = 1.0 - Omega0
HubbleParam = 0.704
num_files = 1
BoxSize = 1.0

mass_array = np.array([0, m[npart[0]], m[npart[0]+npart[1]], m[npart[0]+npart[1]+npart[2]],0,0], dtype=np.float32)

header = g3.GadgetHeader(npart, mass_array, time, redshift, BoxSize, Omega0, OmegaLambda, HubbleParam, num_files=num_files)
f = g3.GadgetWriteFile(output_file, npart, {}, header)
f.write_header(f.header)

f.add_file_block('ID  ', ntot * 4, partlen = 4, dtype=np.int32)
f.add_file_block('POS ', ntot * 4 * 3, partlen = 4 * 3)
f.add_file_block('VEL ', ntot * 4 * 3, partlen = 4 * 3)
f.add_file_block('MASS', npart[0] * 4, partlen = 4)
f.add_file_block('U   ', npart[0] * 4, partlen = 4)

f.write_block('ID  ', -1, id)
f.write_block('POS ', -1, pos)
f.write_block('MASS', -1, m[:npart[0]])
f.write_block('VEL ', -1, vel)
f.write_block('U   ', -1, u[:npart[0]])