import g3read as g3
import numpy as np
from shutil import copyfile

"""
Set up a box with `N` particles per side.
"""
def setup_box(N=20, rho=1.5e-8, P=1.0, 
              gamma=5.0/3.0, boxlen=12000.0):
    
    # total number of particles
    Np = N**3

    # allocate empty arrays
    pos = np.zeros((Np, 3), np.float32)
    vel = np.zeros((Np, 3), np.float32)
    m   = np.zeros(Np, np.float32)
    U   = np.zeros(Np, np.float32)

    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):

                # current active particle 
                ipart = i*N**2 + j*N + k

                # define positions 
                pos[ipart,0] = ( i + 0.5 ) * boxlen / N
                pos[ipart,1] = ( j + 0.5 ) * boxlen / N
                pos[ipart,2] = ( k + 0.5 ) * boxlen / N

                # define velocities 
                vel[ipart,0] = 0.0
                vel[ipart,1] = 0.0
                vel[ipart,2] = 0.0

                # define internal energy
                U[ipart]     = P / ( gamma - 1.0 ) / rho

                # define mass 
                m[ipart]     =  (rho * boxlen**3)/Np

    # set up IDs
    id = np.arange(0, Np, dtype=np.int32)+1

    cond = np.where( (pos[:,0] - boxlen/2)**2 + (pos[:,1] - boxlen/2)**2 + (pos[:,2] - boxlen/2)**2 < (boxlen/2)**2) 
    pos = pos[cond]
    vel = vel[cond]
    U = U[cond]
    m = m[cond]
    id = id[cond] 

    return pos, vel, id, m, U

def write_output(N, pos, vel, id, m, U, output_file):

    # copy reference file to output
    copyfile(r"C:\Users\johan\Python\Hydrodynamics\Tutorial01\ic_header", output_file)

    # read reference file
    
    k = len(pos)
    
    f = g3.GadgetFile(output_file, is_snap=False) #is_snap=False skip the POS block
    
    f.header.npart = [0, k, 0, 0, 0, 0]
    f.header.npartTotal = [0, k, 0, 0, 0, 0]
    f.header.time = 0.0
    f.header.boxsize = 1.0

    f.write_header(f.header,filename=output_file)


    
    f.add_file_block('POS ', (k)*4*3  , partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
    f.add_file_block('VEL ', (k)*4*3, partlen=4*3) #add a block of N^3*4*3 bytes, and each particle takes  4*3 bytes (3 floats)
    f.add_file_block('ID  ', (k)*4, partlen=4, dtype=np.int32) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 int)
    f.add_file_block('MASS', (k)*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)
    #f.add_file_block('U   ', (N**3)*4, partlen=4) #add a block of N^3*4 bytes, and each particle takes  4 bytes (1 float)

    # print(N**3, 'len pos',len(pos), 'len vel',len(vel), 'len IDs',len(id))
    f.write_block("POS ", 1, pos, filename=output_file)
    f.write_block("VEL ", 1, vel, filename=output_file)
    f.write_block("ID  ", 1, id,  filename=output_file)
    f.write_block("MASS", 1, m,   filename=output_file)
    #f.write_block("U   ", 1, U,   filename=output_file)


def main(N, output_file, boxlen):

    # set up the box
    pos, vel, id, m, U = setup_box(N, boxlen = boxlen)

    # write the IC files
    write_output(N, pos, vel, id, m, U, output_file)
    return pos, m


def test(output_file):
    f = g3.GadgetFile(output_file)
    print(output_file,'has','part array of',f.header.npart)

    #readed_ids = f.read_new('ID  ',0)
    #readed_poses = f.read_new('POS ',0)
    #readed_masses = f.read_new('MASS',0)
    
    #print(output_file,'has','ID array of', readed_ids)
    #print(output_file,'has','pos array of', readed_poses)
    #print(output_file,'has','mass array of', readed_masses)
    
    
N = 20
output_file = "box.ic"
pos, mass = main(N, output_file, 12000)
test(output_file)