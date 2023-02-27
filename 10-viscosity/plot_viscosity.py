import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,220,221)

for i in range(len(frames)):
    filename = f'viscosity/visc_{int(frames[i]):03}'
    ptype = 0
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    u = g3.read_new(filename, "U   ", ptype)
    id = g3.read_new(filename, "ID  ", ptype)
    vel = g3.read_new(filename, "VEL ", ptype)

    fig= plt.figure(figsize = (5,9), dpi = 150)
    
    plt.scatter(pos[::4,0], pos[::4,1], marker = '.', c = u[::4], s = 0.1, cmap = 'turbo', vmin = 1, vmax = 11)
    cbar = plt.colorbar()
    cbar.set_label(r'internal energy')
    
    plt.xlabel('x')
    plt.ylabel('y')
    
    
    
    plt.title(f'time = {time:.3f}')
    
    plt.savefig(f'visc_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()
