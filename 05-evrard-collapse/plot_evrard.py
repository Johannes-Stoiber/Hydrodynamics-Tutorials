import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

gamma = 5/3
R_star = 1.0
M_star = 1.5815 
G = 1.0

rho_char = 3*M_star/(4*np.pi*R_star**3)
e_char = G*M_star/R_star
p_char = rho_char*e_char


for i in range(0,10):
    filename = f'snap_00{i}'
    ptype = 0
    time = i*0.2/0.7952

    pos = g3.read_new(filename,"POS ",ptype)
    density = g3.read_new(filename, "RHO ", ptype)/rho_char
    np.random.seed(2)
    np.random.shuffle(pos)
    np.random.seed(2)
    np.random.shuffle(density)

    fig= plt.figure(figsize = (9,7))

    plt.scatter(pos[:,0], pos[:,1], marker = '.', c = density, s = 1, cmap = 'viridis', vmin = 20, vmax = 40)
    cbar = plt.colorbar()
    cbar.set_label(r'$\rho/\rho_{char}$')

    plt.xlabel('x [R_star]')
    plt.ylabel('y [R_star]')
    plt.xlim(-1.,1.)
    plt.ylim(-1.,1.)
    
    plt.text(0.4,0.9,f'time/t_char = {time:.3f}')
    
    plt.savefig(f'evrard_00{i}.png')
    plt.show()