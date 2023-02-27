import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

#frames = np.linspace(0,25,26)
frames = np.linspace(0,25,26)
times = frames*0.2

gamma = 5/3
P = 1.0
rho=0.125


for i in range(len(frames)):
    pos = g3.read_new(f'snap_{int(frames[i]):03}','POS ', 0)
    vel = g3.read_new(f'snap_{int(frames[i]):03}','VEL ', 0)
    density = g3.read_new(f'snap_{int(frames[i]):03}', 'RHO ', 0)
    u = g3.read_new(f'snap_{int(frames[i]):03}', 'U   ', 0)
    pressure = u*(gamma-1)*density
    
    fig = plt.figure( figsize = (10,1), dpi = 150 )
    
    #plt.scatter(pos[::4,0], vel[::4,0], s = 0.1, color = 'm', label = 'v_x' )
    plt.scatter(pos[::4,0], density[::4], s = 0.1, color = 'c', label = 'rho' )
    #plt.scatter(pos[::4,0], pressure[::4], s = 0.1, color = 'b', label = 'P')
    
    plt.legend(loc = 'upper right')
    #plt.text(,, f'{times[i]:.2f}')
    
    #plt.ylim(0.04,0.055)
    
    
    print(len(pos))

    #plt.savefig(f'profiles_{int(frames[i])}.png')
    plt.show()

