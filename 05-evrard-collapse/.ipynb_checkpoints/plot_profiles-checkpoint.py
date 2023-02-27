import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,9,10)
#frames = np.array([10])
color = ['r', 'b', 'g', 'cyan']
times = frames*0.2/0.7952

gamma = 5/3
R_star = 1.0
M_star = 1.5815 
G = 1.0

rho_char = 3*M_star/(4*np.pi*R_star**3)
e_char = G*M_star/R_star
p_char = rho_char*e_char

for i in range(len(frames)):
    pos = g3.read_new(f'snap_00{int(frames[i])}','POS ', 0)
    density = g3.read_new(f'snap_00{int(frames[i])}', 'RHO ', 0)/rho_char
    radii = np.sqrt( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 )/R_star
    u = g3.read_new(f'snap_00{int(frames[i])}', 'U   ', 0)/e_char
    pressure = u*(gamma-1)*density/10
    
    fig = plt.figure( figsize = (5,3), dpi = 150 )
    plt.title( f'time/t_char = {times[i]:.3f}' )
    
    plt.scatter(radii, density, s = 0.1, color = 'm', label = 'density' )
    plt.scatter(radii, pressure, s = 0.5, color = 'c', label = 'pressure' )
    
    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'Radius $r/R_{star}$')
    plt.ylabel(r'$\rho/\rho_{char}$ and $p/p_{char}/10$')

    plt.xticks([0.001, 0.01, 0.1, 1.])
    plt.yticks([0.1, 1.0, 10., 100., 1000.])

    plt.ylim(0.05, 3000)
    plt.xlim(0.001,2)

    plt.legend()
    
    plt.savefig(f'profiles_{int(frames[i])}.png')
    plt.show()