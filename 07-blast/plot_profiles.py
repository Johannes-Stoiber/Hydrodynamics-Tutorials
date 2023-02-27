import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(1,30,30)

gamma = 5./3.
rho = 0.125
Pback = 1e-20
Pblast = 1
for i in range(len(frames)):
    filename = f'blast1/blast_{int(frames[i]):03}'
    ptype = 0
    
    file = g3.GadgetFile(filename)
    time = file.header.time

    pos = g3.read_new(filename,"POS ",ptype)
    density = g3.read_new(filename, "RHO ", ptype)
    vel = g3.read_new(filename, "VEL ", ptype)
    radius = np.sqrt( (pos[:,0]-0.5)**2 + (pos[:,1]-0.5)**2 + (pos[:,2]-0.5)**2 )
    abs_v = np.sqrt( vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2 )
    u = g3.read_new(filename, 'U   ', ptype)
    pressure = u*(gamma-1)*density
    
    ffig = plt.figure(figsize = (11,3), dpi = 200)
    grid = plt.GridSpec(3,11, wspace = 0.0, hspace = 0.0)
    ax1 = plt.subplot(grid[0:3,0:3])
    ax2 = plt.subplot(grid[0:3,4:7])
    ax3 = plt.subplot(grid[0:3,8:11])

    ax1.set_xlabel('radius')
    ax2.set_xlabel('radius')
    ax3.set_xlabel('radius')

    ax1.set_ylabel('pressure')
    ax2.set_ylabel('density')
    ax3.set_ylabel('velocity')

    ax1.scatter(radius[:], pressure[:], marker = '.', color = 'b', s = 0.1)
    ax2.scatter(radius[:], density[:], marker = '.', color = 'm', s = 0.1)
    ax3.scatter(radius[:], abs_v[:], marker = '.', color = 'r', s = 0.1)
    
    
    R = 0.3160358*time**(2/5) 
    ax1.axvline(x = R, color = 'gray', linestyle = '--')
    ax2.axvline(x = R, color = 'gray', linestyle = '--')
    ax3.axvline(x = R, color = 'gray', linestyle = '--')
    
    P = 0.00114911*time**(-6/5)
    ax1.axhline(y = P, color = 'gray', linestyle = '--') 
    
    #rho0 = 0.125
    ax2.axhline(y = rho, color = 'gray', linestyle = '--')
    
    V = 0.10180965*time**(-3/5)
    ax3.axhline(y = V, color = 'gray', linestyle = '--') 
    
    ax1.set_ylim(1e-20,1) #pressure
    ax1.set_yscale('log') #for pressure
    ax2.set_ylim(0,0.4) #rho
    ax3.set_ylim(0,0.5) #velocity
    
    ax1.text(0.5,0.035,f'time= {time:.3f}') #pressure
    ax2.text(0.5,0.37,f'time= {time:.3f}') #rho
    ax3.text(0.5,0.455,f'time= {time:.3f}') #velocity
    
    plt.savefig(f'blast_profiles_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()