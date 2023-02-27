import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

#frames = np.linspace(0,10,11)
frames = np.linspace(0,30,31)
times = frames*0.2

gamma=5./3.
rho = 1.0
P = 1.0

for i in range(len(frames)):
    print(i)
    pos = g3.read_new(f'snaps/mhd_{int(frames[i]):03}','POS ', 0)
    vel = g3.read_new(f'snaps/mhd_{int(frames[i]):03}','VEL ', 0)
    density = g3.read_new(f'snaps/mhd_{int(frames[i]):03}', 'RHO ', 0)
    u = g3.read_new(f'snaps/mhd_{int(frames[i]):03}', 'U   ', 0)
    b = g3.read_new(f'snaps/mhd_{int(frames[i]):03}', 'BFLD', 0)
    divb = g3.read_new(f'snaps/mhd_{int(frames[i]):03}', 'DIVB', 0)
    pressure = u*(gamma-1)*density
    
    contact_disc = np.where(divb == np.max(divb))
    x_contact_disc = pos[contact_disc[0],0]
    
    
    fig = plt.figure(figsize = (10,6), dpi = 300)
    grid = plt.GridSpec(3,3, wspace = 0.4, hspace = 0.3)
    ax1 = plt.subplot(grid[0:1,0:1])
    ax2 = plt.subplot(grid[1:2,0:1])
    ax3 = plt.subplot(grid[2:3,0:1])
    ax4 = plt.subplot(grid[0:1,1:2])
    ax5 = plt.subplot(grid[1:2,1:2])
    ax6 = plt.subplot(grid[2:3,1:2])
    ax7 = plt.subplot(grid[0:1,2:3])
    ax8 = plt.subplot(grid[1:2,2:3])
    ax9 = plt.subplot(grid[2:3,2:3])
    
    ax1.scatter(pos[::4,0], density[::4], s = 0.1, color = 'c', label = 'rho' )
    ax1.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    ax2.scatter(pos[::4,0], vel[::4,0], s = 0.1, color = 'red', label = 'vel_x' )
    ax2.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    ax3.scatter(pos[::4,0], divb[::4], s = 0.1, color = 'blue', label = 'div(B)' )
    ax3.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    
    
    ax4.scatter(pos[::4,0], pressure[::4], s = 0.1, color = 'c', label = 'pressure' )
    ax4.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    ax5.scatter(pos[::4,0], vel[::4,1], s = 0.1, color = 'red', label = 'vel_y' )
    ax5.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    ax6.scatter(pos[::4,0], b[::4,1], s = 0.1, color = 'blue', label = 'B_y' )
    ax6.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    
    ax7.scatter(pos[::4,0], u[::4], s = 0.1, color = 'c', label = 'internal energy' )
    ax7.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    ax8.scatter(pos[::4,0], vel[::4,2], s = 0.1, color = 'red', label = 'vel_z' )
    ax8.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    ax9.scatter(pos[::4,0], b[::4,2], s = 0.1, color = 'blue', label = 'B_z' )
    ax9.axvline(x = x_contact_disc, linestyle = '--', color = 'gray', linewidth = 0.5)
    
    
    ax1.set_ylabel('density')
    ax2.set_ylabel(r'$vel_x$')
    ax3.set_ylabel(r'$div(B)$')
    
    ax4.set_ylabel('pressure')
    ax5.set_ylabel(r'$vel_y$')
    ax6.set_ylabel(r'$B_y$')
    
    ax7.set_ylabel('internal energy')
    ax8.set_ylabel(r'$vel_z$')
    ax9.set_ylabel(r'$B_z$')
    
    ax1.set_xlim(30,70)
    ax2.set_xlim(30,70)
    ax3.set_xlim(30,70)
    ax4.set_xlim(30,70)
    ax5.set_xlim(30,70)
    ax6.set_xlim(30,70)
    ax7.set_xlim(30,70)
    ax8.set_xlim(30,70)
    ax9.set_xlim(30,70)
    
    ax1.set_ylim(-0.2,1.2)
    ax2.set_ylim(-0.4,0.8)
    ax3.set_ylim(-1.2,1.2)
    ax4.set_ylim(-0.1,1.2)
    ax5.set_ylim(-2,0.5)
    ax6.set_ylim(-1.2,1.2)
    ax7.set_ylim(0,4)
    ax8.set_ylim(-1.2,1.2)
    ax9.set_ylim(-1.2,1.2)
    
    ax4.set_title(f'time = {times[i]:.2f}')
    
    plt.savefig(f'mhd_{int(frames[i]):03}.png', format = 'png', bbox_inches = 'tight')
    plt.show()