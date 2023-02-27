import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,107,108)

for i in range(len(frames)):
    filename = f'output/snap_{int(frames[i]):03}'
    ptype = 0
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    u = g3.read_new(filename, "U   ", ptype)
    id = g3.read_new(filename, "ID  ", ptype)

    fig= plt.figure(figsize = (12,2), dpi = 250)
    
    boundary = np.where(id[:] == 0)
    free = np.where(id[:] != 0)
    
    plt.scatter(pos[free][:,0], pos[free][:,1], marker = '.', c = u[free], s = 0.1, cmap = 'turbo', vmin = 8, vmax = 12)
    cbar = plt.colorbar()
    cbar.set_label(r'internal energy')
    
    plt.scatter(pos[boundary][::2,0], pos[boundary][::2,1], s = 0.5, color = 'gray', label = 'boundary particles')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(10,40)
    
    plt.title(f'nozzle at time = {time:.3f}')
    
    plt.savefig(f'nozzle_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()