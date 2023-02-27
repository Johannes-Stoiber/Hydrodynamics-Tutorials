import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

filename = r'C:\Users\johan\Python\Hydrodynamics\Tutorial01\snap_050'
ptype = 1

pos = g3.read_new(filename,"POS ",ptype)

fig = plt.figure(figsize = (7,7))
plt.plot(pos[:,0], pos[:,1], '.')
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
plt.xlim(0,12000)
plt.ylim(0,12000)
plt.savefig('snap_050.jpg')
plt.show()