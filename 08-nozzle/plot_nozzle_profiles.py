import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,107,108)

gamma = 5./3.
rho = 0.125
P = 1.0
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
    
    fig= plt.figure(figsize = (12,2), dpi = 250)

    #plt.scatter(pos[free][:,0], pressure[free][:], marker = '.', color = 'm', s = 0.1)
    #plt.scatter(pos[wall][:,0], pressure[wall][:], marker = '.', color = 'gray', s = 0.5)
    plt.scatter(pos[free][:,0], vel[free][:,0], marker = '.', color = 'm', s = 0.1)
    plt.scatter(pos[wall][:,0], vel[wall][:,0], marker = '.', color = 'gray', s = 0.5)
    
    plt.scatter(pos[free][:,0], -sound_speed[free][:], marker = '.', color = 'r', s = 0.1)
    plt.scatter(pos[wall][:,0], -sound_speed[wall][:], marker = '.', color = 'gray', s = 0.5)
    
    
    plt.axvline(x = min(pos[wall][:,0]), color = 'gray', linestyle = '--')
    plt.axvline(x = min(pos[wall][:,0])+0.5, color = 'gray', linestyle = '--')
    

    plt.xlabel('x')
    plt.ylabel(r'$c_s$ and $v_x$')
    
    plt.xlim(10,40)
    
    #plt.ylim(0,0.5) #pressure
    plt.ylim(-5,0.5) #velocity
    
    plt.title(f'velocity in x-direction at time = {time:.3f}')
    
    plt.savefig(f'nozzle_velocity_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()