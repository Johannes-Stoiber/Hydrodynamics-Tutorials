import numpy as np
import matplotlib.pyplot as plt

#N = np.array([16,32,48,64,80,96,112,128,144,160,176,192,208,224,240,256])
n = 16
N = []
while n < 257:
    N.append(n)
    n += 8

kernel = ["wend_c2",  "wend_c6", "cubic", "quintic"]

density = np.zeros((len(kernel), len(N))) 

for i in range(len(kernel)):
    for j in range(len(N)):
        text = open("density_{}_{}.txt".format(kernel[i], N[j]))
        zeilen = text.readlines()[0]
        zi = zeilen.split()
        n = zi[0]
        ratio = zi[4]
        density[i,j] = ratio

log_density = np.log(density)


fig = plt.figure(figsize = (9,6))
for i in range(len(kernel)):
    plt.plot(N, log_density[i], '.', label = kernel[i])

plt.xscale('log') 
plt.xlabel('N')
plt.ylabel('log(rhosim/rhoini)')
plt.legend(loc = 0)
plt.grid()
plt.savefig('densities.png')
plt.show()
    
