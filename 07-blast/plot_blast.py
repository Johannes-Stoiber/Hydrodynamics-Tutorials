import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,30,31)

for i in range(len(frames)):
    filename = f'blast2/blast2_{int(frames[i]):03}'
    ptype = 0
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    density = g3.read_new(filename, "RHO ", ptype)

    fig= plt.figure(figsize = (9,7.2))
    
    cond = np.where((pos[:,2] < 0.55) & (pos[:,2] > 0.45))
    pos2 = np.copy(pos[cond])

    plt.scatter(pos2[:,0], pos2[:,1], marker = '.', c = density[cond], s = 1, cmap = 'viridis', vmin = 0.048, vmax = 0.32 )
    cbar = plt.colorbar()
    cbar.set_label(r'density')

    plt.xlabel('x')
    plt.ylabel('y')
    #plt.xlim(-1.,1.)
    #plt.ylim(-1.,1.)
    
    #plt.text(0.4,0.9,f'time= {time:.3f}')
    plt.title(f'blast wave at time = {time:.3f}')
    
    plt.savefig(f'blast_new_{int(frames[i]):03}.png', format = 'png')
    plt.show()