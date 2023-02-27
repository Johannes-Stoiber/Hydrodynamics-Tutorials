# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 17:22:03 2022

@author: user
"""

# -*- coding: utf-8 -*-
"""
Hydrodynamics Tutorial 5
Visualsation For Eint, Epot, Etot

Created on Mon Nov 21 17:30:59 2022

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

file = "energy.txt"

# read the energies
data = np.genfromtxt(file, comments="%", dtype=float)
e_char = 1.5815
t_char = 0.7952

# plot all things
time = data[:,0] / t_char
Eint = data[:,1] / e_char
Epot = data[:,2] / e_char
Ekin = data[:,3] / e_char
Etot = Eint + Epot + Ekin
ratio = Eint/Epot

fig = plt.figure( figsize = (9,6))

plt.scatter(time, Eint, color = 'r', label="Eint", s=0.5)
plt.scatter(time, Ekin, color = 'b', label="Ekin", s=0.5)
plt.scatter(time, Etot, color = 'k', label="Etot", s=0.5)
plt.scatter(time, Epot, color = 'g', label="Epot", s=0.5)
plt.scatter(time[:], ratio[:], color = 'm', label = 'virial ratio', s = 0.5)

plt.axhline(y=0, linestyle = '--', color = 'gray')
plt.axhline(y=-0.5, linestyle = '--', color = 'gray')

plt.legend(fontsize = 12)
plt.xlabel("$t/t_{char}$", fontsize = 12)
plt.ylabel("$E/e_{char}$", fontsize = 12)

plt.yticks(fontsize = 12)
plt.xticks(fontsize = 12)

plt.savefig("Energyplot.png")
plt.show()

