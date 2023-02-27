import g3read as g3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col

frames = np.linspace(0,107,108)

gamma = 5./3.
rho = 0.125
Pback = 1e-20
Pblast = 1
for i in range(len(frames)):
    filename = f'output/snap_{int(frames[i]):03}'
    ptype = 0
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    density = g3.read_new(filename, "RHO ", ptype)
    vel = g3.read_new(filename, "VEL ", ptype)
    u = g3.read_new(filename, 'U   ', ptype)
    id = g3.read_new(filename, 'ID  ', ptype)
    pressure = u*(gamma-1)*density
    sound_speed = np.sqrt(gamma*pressure/density)
    
    boundary = np.where(id[:] == 0)
    wall = np.where((id[:] == 0) & (vel[:,0] == -1))
    free = np.where(id[:] != 0)
    
    fig = plt.figure(figsize = (12,6), dpi = 200)
    grid = plt.GridSpec(8,1, wspace = 0.0, hspace = 0.0)
    ax1 = plt.subplot(grid[0:2,0:1])
    ax2 = plt.subplot(grid[3:5,0:1])
    ax3 = plt.subplot(grid[6:8,0:1])

    ax1.set_xlabel('x')
    ax2.set_xlabel('x')
    ax3.set_xlabel('x')

    ax1.set_ylabel('z')
    ax2.set_ylabel('P')
    ax3.set_ylabel(r'$c_s$ and $v_x$')
    
    ax1.scatter(pos[free][:,0], pos[free][:,1], marker = '.', c = u[free], s = 0.1, cmap = 'turbo', vmin = 8, vmax = 12)
    cbar = plt.colorbar(cm.ScalarMappable(norm = col.Normalize(vmin = 8, vmax = 12),cmap= 'turbo'), ax=ax1, location = 'top')
    #cbar = plt.colorbar(cax = ax1)
    cbar.set_label(r'internal energy')
    ax1.scatter(pos[boundary][::2,0], pos[boundary][::2,1], s = 0.5, color = 'gray', label = 'boundary particles')
    
    ax2.scatter(pos[free][:,0], pressure[free][:], marker = '.', color = 'm', s = 0.1)
    ax2.scatter(pos[wall][:,0], pressure[wall][:], marker = '.', color = 'gray', s = 0.5)
    ax2.axvline(x = min(pos[wall][:,0]), color = 'gray', linestyle = '--')
    ax2.axvline(x = min(pos[wall][:,0])+0.5, color = 'gray', linestyle = '--')
    
    ax3.scatter(pos[free][:,0], vel[free][:,0], marker = '.', color = 'm', s = 0.1)
    ax3.scatter(pos[wall][:,0], vel[wall][:,0], marker = '.', color = 'gray', s = 0.5)
    ax3.scatter(pos[free][:,0], -sound_speed[free][:], marker = '.', color = 'r', s = 0.1)
    ax3.scatter(pos[wall][:,0], -sound_speed[wall][:], marker = '.', color = 'gray', s = 0.5)
    ax3.axvline(x = min(pos[wall][:,0]), color = 'gray', linestyle = '--')
    ax3.axvline(x = min(pos[wall][:,0])+0.5, color = 'gray', linestyle = '--')
    
    ax1.set_ylim(0,1) #flow
    ax2.set_ylim(0,0.5) #pressure
    ax3.set_ylim(-5,0.5) #velocity
    
    ax1.set_xlim(10,40)
    ax2.set_xlim(10,40)
    ax3.set_xlim(10,40)
    
    #ax1.set_title(f'nozzle at time = {time:.3f}')
    ax2.set_title(f'pressure at time = {time:.3f}')
    ax3.set_title(f'velocity in x-direction at time = {time:.3f}')
    
    plt.tight_layout()
    
    plt.savefig(f'nozzle_all_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()