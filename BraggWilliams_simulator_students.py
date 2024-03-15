#!/usr/bin/env python
from __future__ import division
#%matplotlib inline
import matplotlib
from numba import njit
# Simulating the Bragg-Williams model
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

# Interaction parameters
waa = -1.
wbb = -1.
wab = 1.

# System size
N = 500  # sites in x and y-directions
M = N*N # total number of sites



#    ''' Simulating the Bragg-Williams model a la Martin Trulsson'''  
# Code re-adated from https://github.com/rajeshrinet/compPhy and the Ising model
#    ## monte carlo moves
@njit
def mcmove(config, N, beta):
    ''' This is to execute the monte carlo moves using 
    Metropolis algorithm such that detailed
    balance condition is satisified'''

    for i in range(N):
        for j in range(N):            
                x1  = np.random.randint(0, N)
                y1  = np.random.randint(0, N)
                x2 = np.random.randint(0, N)
                y2 = np.random.randint(0, N)
                s1 = config[x1, y1]
                s2 = config[x2, y2]
                if s1!=s2: # only if 1 and 2 is different otherwise does not change state
                    # number of a neighbours place 1
                    nb1 = config[(x1+1)%N,y1] + config[x1,(y1+1)%N] + config[(x1-1)%N,y1] + config[x1,(y1-1)%N]
                    # number of a neighbours place 2 -> What does %N do? 
                    nb2 = config[(x2+1)%N,y2] + config[x2,(y2+1)%N] + config[(x2-1)%N,y2] + config[x2,(y2-1)%N]
                    na1 = 4.-nb1
                    na2 = 4.-nb2
                    if s1==1: # inital 1:b 2:a center # final 1:a 2:b
                        de = -(waa+wbb-2.*wab)*(nb1-nb2)
                    else: # a center
                        de = (waa+wbb-2.*wab)*(nb1-nb2)
                    if de < 0: # energy decreases -> always accept	
                        if s1==1: 
                            config[x1,y1]=0
                            config[x2,y2]=1
                        else:
                            config[x1,y1]=1
                            config[x2,y2]=0
                    elif np.random.random() < np.exp(-de*beta): 
                        if s1==1:
                            config[x1,y1]=0
                            config[x2,y2]=1
                        else:
                            config[x1,y1]=1
                            config[x2,y2]=0
    return config
    
def simulate():   
    ''' This module simulates the Bragg-Williams model'''
    # Initialse the lattice (random configuration)
    config = np.zeros((N,N))   # Matrix N*N
    xmix = 0.5
    # In the while loop, a fraction of the lattice sites are assigned to particle B (value set to 1)
    i = 0
    while i < M*xmix:
        x  = np.random.randint(0, N)
        y  = np.random.randint(0, N)
        if config[x,y]==0:
          config[x,y]=1
          i=i+1

    # Plot data
    f = plt.figure(figsize=(15, 15), dpi=80);    
    plt.suptitle("Temperature="+str(temp))
    # Plot the first random configuartion
    configPlot(f, config, 0, N, 1);
    msrmnt = 50001
    intU = np.empty(0)
    for i in range(msrmnt):
        mcmove(config, N, 1.0/temp)
        if i == 1:       configPlot(f, config, i, N, 2);
        if i == 50:       configPlot(f, config, i, N, 3);
        if i == 1000:      configPlot(f, config, i, N, 4);
        if i == 20000:     configPlot(f, config, i, N, 5);
        if i == 50000:     configPlot(f, config, i, N, 6);
        # Calculate the total energy
        if i in np.arange(200,500,1):
            u = 0
            for i in range(N):
                for j in range(N): 
                    # Find an expression for the total energy
                    u=u+0.
            intU = np.append(intU,u)         
    return intU.mean()/(M),intU.std()/float(M)
                 
                    
def configPlot(f, config, i, N, n_):
    ''' This modules plts the configuration once passed to it along with time etc '''
    X, Y = np.meshgrid(range(N), range(N))
    sp =  f.add_subplot(3, 3, n_ )  
    plt.setp(sp.get_yticklabels(), visible=False)
    plt.setp(sp.get_xticklabels(), visible=False)      
    plt.imshow(config, interpolation='none',cmap=plt.cm.RdBu)
    plt.title('Time=%d'%i)     

intUvsT = np.empty(0)
intUerrvsT = np.empty(0)
intU2vsT = np.empty(0)
tempvec = np.arange(6.,6.5,0.5)
for temp in tempvec:
    # rm = BraggWilliams()
    intU_mean, intU_std = simulate()
    intUvsT = np.append(intUvsT, intU_mean)
    intUerrvsT = np.append(intUerrvsT, intU_std)
    
plt.show()
plt.errorbar(tempvec, intUvsT/tempvec, intUerrvsT/tempvec, lw=0, marker='o',color='black')
plt.title('Internal Energy')
plt.ylabel('$U$  / $kT$')
plt.xlabel('$kT$')
plt.show()
