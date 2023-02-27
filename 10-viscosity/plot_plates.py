import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,3,4)

for i in range(len(frames)):
    filename = f'plates/plates_{int(frames[i]):03}'
    ptype = 0
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    u = g3.read_new(filename, "U   ", ptype)
    id = g3.read_new(filename, "ID  ", ptype)
    vel = g3.read_new(filename, "VEL ", ptype)
    
    #vel0 = np.sqrt( vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2 )

    fig= plt.figure(figsize = (5,4.5), dpi = 250)
    
    boundary = np.where(id[:] == 0)
    free = np.where(id[:] != 0)
    
    plt.scatter(pos[::4,0], pos[::4,1], marker = '.', c = vel[::4,1], s = 0.1, cmap = 'turbo')
    cbar = plt.colorbar()
    cbar.set_label(r'velocity')
    
    plt.scatter(pos[::4,0], vel[::4,1]/2+0.5, s = 0.2, color = 'black')
    
    plt.xlabel('x')
    plt.ylabel('y')
    
    plt.title(f'time = {time:.3f}')
    
    #plt.savefig(f'plates_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()