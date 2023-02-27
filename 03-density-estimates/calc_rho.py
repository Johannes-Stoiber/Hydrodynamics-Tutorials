import g3read1 as g3
import numpy as np
import sys

N = sys.argv[1]
kernel = "wend_c6"
filename = "snap_000"

text = open("density_{}_{}.txt".format(kernel,N), "w")
#f.write("n \t rhoini \t rhosim \t shovar \t ratio \n")

f = g3.GadgetFile(filename)
boxsize = f.header.BoxSize
npart = f.header.npart[0]

pos = g3.read_new(filename, "POS ", 0)
rho = g3.read_new(filename, "RHO ", 0)
mass = g3.read_new(filename, "MASS", 0)

rhoini = np.sum(mass)/boxsize**3
rhosim = np.sum(rho)/npart
rhovar = np.sum((rho-rhosim)**2)/npart
ratio = rhosim/rhoini
    
text.write(str(N)+"\t"+str(rhoini)+"\t"+str(rhosim)+"\t"+str(rhovar)+"\t"+str(ratio)+"\n")
text.close() 
