# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 14:07:07 2019

@author: YYX
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import random
#parameter  
lambd = 1  # step size is 1
random.seed(151)

# (1)
plt.figure(1)
x = 0
y = 0
n = 100
# calculate Distance for a given steps
for i in range(n):
    theta = random.random()*2*np.pi
    dx = lambd*np.cos(theta)
    dy = lambd*np.sin(theta)
    if i != n-1:
        plt.arrow(x, y, dx, dy,length_includes_head=True,
                  head_width=0.04, head_length=0.04)
    # the last arrow is marked red    
    if i == n-1:
        plt.arrow(x, y, dx, dy,length_includes_head=True,
                  head_width=0.04, head_length=0.04, color = 'r') 
    x += dx
    y += dy
    
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size =10)
plt.xlabel('x')
plt.ylabel('y') 
plt.axis([-15,15,-15,15])
plt.savefig('C:/Users/YYX/Desktop/ANU/PHYS3103/comp_project/example_walk.pdf')
    


# (2)
D = []
N = np.logspace(1,6,1000)
N = N.astype(int)
plt.figure(2)

for n in N:
    DD = 0
    for j in range(10):   
        # calculate Distance for a given steps
        x = 0
        y = 0
        for i in range(n):
            theta = random.random()*2*np.pi
            x += lambd*np.cos(theta)
            y += lambd*np.sin(theta)
        DD = np.sqrt((DD**2*j+x**2+y**2)/(j+1))   
    D.append(DD)
        
    
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size =10)
plt.xlabel('Steps')
plt.ylabel('Distance')
plt.loglog(N,D, color='navy')
plt.savefig("C:/Users/YYX/Desktop/ANU/PHYS3103/comp_project/Distance_steps.pdf")




