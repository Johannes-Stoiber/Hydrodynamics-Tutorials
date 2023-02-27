import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,0,1)

for i in range(len(frames)):
    #filename = f'galaxy_{int(frames[i]):03}'
    filename = 'galaxy.ic'
    ptype = [0,2,3]
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    ptype = [0,2,3]

    pos = g3.read_new(filename,"POS ",ptype)
    vel = g3.read_new(filename,"VEL ",ptype)

    fig = plt.figure(figsize = (7.5,3), dpi = 150)
    grid = plt.GridSpec(1,2, wspace = 0.4, hspace = 0.2)
    ax1 = plt.subplot(grid[0:1,0:1])
    ax2 = plt.subplot(grid[0:1,1:2])

    #ax1.scatter(pos[:,1], pos[:,2], marker = '.', c = vel[:,0], s = 1, cmap = 'jet_r')
    ax1.scatter(pos[:,0], pos[:,1], marker = '.', c = vel[:,2], s = 0.5, cmap = 'jet_r', vmin = -200, vmax = 200)
    ax2.scatter(pos[:,1], pos[:,2], marker = '.', c = vel[:,0], s = 0.5, cmap = 'jet_r', vmin = -200, vmax = 200)
    #cbar = plt.colorbar()
    #ax2.set_label(r'$v_{los}$')

    ax1.set_xlabel('x [kpc]')
    ax1.set_ylabel('y [kpc]')
    ax1.set_xlim(-40.,40.)
    ax1.set_ylim(-40.,40.)
    
    ax2.set_xlabel('x [kpc]')
    ax2.set_ylabel('z [kpc]')
    ax2.set_xlim(-40.,40.)
    ax2.set_ylim(-40.,40.)
    
    ax1.set_title(f'face on at time = {time:.3f}')
    ax2.set_title(f'edge on at time = {time:.3f}')
    
    #plt.savefig(f'galaxy_{int(frames[i]):03}.png', format = 'png')
    plt.show()