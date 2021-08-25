# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 21:36:02 2018

@author: yyx
"""

# Program to numerically solve the Friedmann equation.
#
# Note - in our derivation of the Friedmann equation we assumed that
# the sphere of the universe was growing (i.e. "a dot" was positive).
# There is an equally valid solution with "a dot" negative, and this
# program is set up to realise when we need to flip to that different
# (contracting universe) solution.

# Program uses Euler's method and the differential form of the equation
# for simplicity sake. A more accurate solution can be obtained using 
# a higher order numerical approximation and the integral form of 
# Friedmann's equation. This one gives perfectly acceptable answers, however.

# Program assumes a matter dominated universe.

# Paul Francis, 2/10/14

import pylab # This program requires the matplotlib python library for the plots.
import math


# Set up constants
G=6.67e-11
c=3.0e8
pc = 3.08e16
pi = math.pi
yearsec = 365.24*24*3600.0

# Pick cosmological model: edit these values as required
H0=70.0 # in km/s/Mpc
Omega_m = 0.3 # Fraction of the critical density in matter
Omega_l = 0.0 # Fraction of the critical density in dark energy
Omega_r = 0.0 # Fraction of the critical density in radiation
# Time range to plot (in years)
trange = 3.0e10

#Compute derived quantities
H0 = H0*1000/(1.0e6*pc) # Convert to SI units
rho0 = 3.0*H0*H0/(8.0*pi*G)
k = (Omega_m + Omega_l-1.0)*H0*H0

#Iterate forward and back

tstep = 1.0e8*yearsec 

# Start going forward from present
a = 1.0

#set up arrays to store times and scale factors
t2arr = []
a2arr = []

# Start loop
t = 0.0
expanding = 1
while (t/yearsec < trange) & (a > 0.0):
    # Work out the densities of matter and Lambda at this time
    rho_m  = Omega_m*rho0/(a*a*a)
    rho_r = Omega_r*rho0/(a*a*a*a)
    rho_l = Omega_l*rho0
    # Work out the rate of change of a(t)
    temp = (8.0/3.0)*pi*G*(rho_m+rho_l+rho_r)-k/(a*a)
    if (temp < 0.0) & expanding == 1: #need different version of Friedmann
        expanding = 0               # equation for expanding and contracting
        temp = -1.0*temp
        
    if expanding == 1:
        adot = a*math.sqrt(temp)
    else: # Equation needs slight modification for contracting universe
        adot = -1.0*a*math.sqrt(temp)
    # work out the new a(t) and t
    a+= adot*tstep
    t += tstep
    
    # Store the new values
    t2arr.append(t/(1.0e9*yearsec))
    a2arr.append(a)
    
# Now go backwards from present    
a = 1.0

#set up arrays to store times and scale factors
t1arr = []
a1arr = []

# Start loop
t = 0.0
expanding  = 1
while (t/yearsec > -1.0*trange) & (a > 0.0):
    # Work out the densities of matter and Lambda at this time
    rho_m  = Omega_m*rho0/(a*a*a)
    rho_r = Omega_r*rho0/(a*a*a*a)
    rho_l = Omega_l*rho0
    # Work out the rate of change of a(t)
    temp = (8.0/3.0)*pi*G*(rho_m+rho_l+rho_r)-k/(a*a)

    if (temp < 0.0) & expanding == 1: #need different version of Friedmann
        expanding = 0               # equation for expanding and contracting
        temp = -1.0*temp
        
    if expanding == 1:
        adot = a*math.sqrt(temp)
    else: # Equation needs slight modification for contracting universe
        adot = -1.0*a*math.sqrt(temp)
        
    # work out the new a(t) and t
    a = a - adot*tstep
    t = t - tstep   
     
    # Store the new values
    t1arr.append(t/(1.0e9*yearsec))
    a1arr.append(a)    
    
# Combine time going forward and back into single arrays

l1 = pylab.size(t1arr)
l2 = pylab.size(t2arr)

tarr = pylab.zeros(l1+l2)
aarr = pylab.zeros(l1+l2)

maxa = 0.0

for i in range(0,l1+l2):
    if i < l1:
        tarr[i] = t1arr[l1-i-1]
        aarr[i] = a1arr[l1-i-1]
    else:
        tarr[i] = t2arr[i-l1]
        aarr[i] = a2arr[i-l1]
    
    # Record maximum a(t) to use in scaling output
    if aarr[i] > maxa:
        maxa = aarr[i]
    
pylab.plot(tarr,aarr,"g-")
pylab.ylim(0.0,maxa)
pylab.xlabel("Time (Gyr)")
pylab.ylabel("a(t)")
pylab.show()