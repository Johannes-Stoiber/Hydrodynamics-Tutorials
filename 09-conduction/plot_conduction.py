import g3read as g3
import numpy as np
import matplotlib.pyplot as plt

frames = np.linspace(0,100,101)
times = frames*0.02
color = ['r','m', 'm', 'm', 'm', 'm', 'm' ,'m', 'm', 'm','b']
opaque = [1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1]

gamma=5./3.
rho = 1.0
P = 1.0

total_entropy = np.zeros(len(frames))

fig = plt.figure( figsize = (10,4), dpi = 200 )

for i in range(int(len(frames)/10)+1):
    j = 10*i
    pos = g3.read_new(f'snaps/snap_{int(frames[j]):03}','POS ', 0)
    u = g3.read_new(f'snaps/snap_{int(frames[j]):03}', 'U   ', 0)
    
    plt.scatter(pos[::4,0], u[::4], s = 0.1, color = color[i], alpha = opaque[i])

    plt.title(f'internal energy distribution for time = {times[j]:.2f}')
    plt.ylabel(r'$u$')
    plt.xlabel(r'$x$')
    
    plt.ylim(0,8)
    #plt.savefig(f'entropy_dist_{int(frames[j]):03}.png', format = 'png', bbox_inches = 'tight')
plt.show()