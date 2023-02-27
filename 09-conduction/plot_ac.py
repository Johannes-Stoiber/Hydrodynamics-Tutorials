import g3read as g3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col

frames = np.linspace(0,100,101)

gamma = 5./3.
rho = 0.125
P= 1
for i in range(len(frames)):
    filename = f'snaps2/ac_{int(frames[i]):03}'
    ptype = 0
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    u = g3.read_new(filename, 'U   ', ptype)
    id = g3.read_new(filename, 'ID  ', ptype)
    
    boundary = np.where(id[:] == 0)
    wall = np.where((id[:] == 0) & (pos[:,0] > 15) & (pos[:,0] < 19.5))
    free = np.where(id[:] != 0)
    
    fig = plt.figure(figsize = (12,4), dpi = 200)
    grid = plt.GridSpec(5,1, wspace = 0.0, hspace = 0.0)
    ax1 = plt.subplot(grid[0:2,0:1])
    ax2 = plt.subplot(grid[3:5,0:1])

    ax1.set_xlabel('x')
    ax2.set_xlabel('x')

    ax1.set_ylabel('y')
    ax2.set_ylabel(r'$u$')
    
    ax1.scatter(pos[free][:,0], pos[free][:,1], marker = '.', c = u[free], s = 0.1, cmap = 'turbo', vmin = 11.4, vmax = 12.7)
    cbar = plt.colorbar(cm.ScalarMappable(norm = col.Normalize(vmin = 11.4, vmax = 12.7),cmap= 'turbo'), ax=ax1, location = 'top')
    #cbar = plt.colorbar(cax = ax1)
    cbar.set_label(r'internal energy')
    ax1.scatter(pos[boundary][::2,0], pos[boundary][::2,1], s = 0.5, color = 'gray', label = 'boundary particles')
    
    ax2.scatter(pos[wall][::4,0], u[wall][::4], marker = '.', color = 'gray', s = 0.5)
    #ax2.scatter(pos[free][::4,0], u[free][::4], marker = '.', c = u[free][::4], cmap = 'turbo', s = 0.1, vmin = 11.6, vmax = 12.5)
    ax2.scatter(pos[free][::4,0], u[free][::4], marker = '.', c = 'b', s = 0.1)
    ax2.axhline(y = 12.0, color = 'm', linestyle = '--')
    
    ax1.set_ylim(0,1) #flow
    ax2.set_ylim(11.4,12.7)  #internal energy
    
    ax1.set_xlim(0,20)
    ax2.set_xlim(0,20)
    
    #ax1.set_title(f'nozzle at time = {time:.3f}')
    ax2.set_title(f'internal energy at time = {time:.3f}')
    
    #plt.savefig(f'ac_all_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()