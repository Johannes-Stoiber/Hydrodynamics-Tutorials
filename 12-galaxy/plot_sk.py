import g3read as g3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

filename = 'snap_010'
file = g3.GadgetFile(filename)
h = file.header
pos = g3.read_new(filename, 'POS ', 0)
mass = g3.read_new(filename, 'MASS', 0)
sfr = g3.read_new(filename, 'SFR ', 0)

radius = np.sqrt(pos[:,0]**2+pos[:,1]**2)
N = 20
radius_tab = 20*10**(np.arange(N)/(N-1)*3-3)
#radius_tab = np.logspace(0.1,1.3,20)

rho_tab = np.arange(N-1, dtype = np.float32)
sfr_tab = np.arange(N-1, dtype = np.float32)

for i in range(N-1):
    cond = np.where((radius >  radius_tab[i]) & (radius < radius_tab[i+1]))
    
    area = np.pi*(radius_tab[i+1]**2-radius_tab[i]**2)
    rho_tab[i] = np.sum( mass[cond]*1e10)/area/1e6
    sfr_tab[i] = np.sum( sfr[cond] )/area
    
fig = plt.figure( figsize = (8,4), dpi = 120 )

plt.plot(rho_tab, sfr_tab, 'b.', label = 'data') 

sort = np.argsort(rho_tab)

rho_tab = rho_tab[sort]
sfr_tab = sfr_tab[sort]

def sk_relation(rho, a, b):
    return a*rho**b

popt, pcov = curve_fit(sk_relation, rho_tab, sfr_tab, p0 = [1,1])

msg = r' $\rho_{sfr} = a\cdot \rho_{gas}^n$'

plt.plot(rho_tab, sk_relation(rho_tab, *popt), 'r-', label = 'fit'+msg)

plt.xlabel('rho_gas')
plt.ylabel('rho_sfr')

plt.xscale('log')

plt.xlim(2,2.2e3)

plt.legend(title = f'a = {popt[0]:3.2}, n = {popt[1]:3.2}'  , loc = 'upper left')
plt.title('Schmidt-Kennicutt Relation in Galaxy Simulation')

plt.show()