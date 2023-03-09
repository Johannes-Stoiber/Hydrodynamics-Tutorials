import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,51,52)
#frames = np.linspace(0,0,1)

for i in range(len(frames)):
    filename = f'snaps_same/same_{int(frames[i]):03}'
    #filename = 'galaxies_impact.ic'
    ptype = [0,2,3]
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    vel = g3.read_new(filename,"VEL ",ptype)
    
    pos_dm = g3.read_new(filename,"POS ", 1)
    
    sort_ind = np.argsort(-pos[:,0])
    pos = pos[sort_ind]
    vel = vel[sort_ind]

    fig= plt.figure(figsize = (25,8), dpi = 250)
    plt.scatter(pos_dm[:,1], pos_dm[:,2], marker = '.',s = .5, c = 'black')
    plt.scatter(pos[:,1], pos[:,2], marker = '.', c = -vel[:,0], s = .5, cmap = 'coolwarm_r', vmin = -200, vmax = 200)
    
    #plt.scatter(pos_dm[:,1], pos_dm[:,0], marker = '.',s = .5, c = 'gray')
    #plt.scatter(pos[:,1], pos[:,0], marker = '.', c = vel[:,2], s = .5, cmap = 'jet_r', vmin = -200, vmax = 200)
    
    cbar = plt.colorbar()
    cbar.set_label(r'$v_{los} [km/s]$')

    #plt.xlim(-100,100)
    plt.ylim(-30,30)
    plt.xlim(-75,75)
    #plt.ylim(-600,600)
    
    plt.xlabel('y [kpc]')
    plt.ylabel('z [kpc]')
    
    plt.title(f'merger at time = {time:.3f}')
    
    plt.savefig(f'galaxies_same/galaxies_same_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    #plt.show()
