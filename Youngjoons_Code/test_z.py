#!/usr/bin/env python2.7

import matplotlib.pyplot as plt
import numpy as np
from atm import *
import libatm as latm


"""**************************************************************
To call the function data_f90 and mesh, we just define following 
tempoaray data. 
**************************************************************"""
Nx=400
Np=500

x0=0.0
xf=45000.0
pA=10.0
p0=1000.0


DEBUG=0
dt=1
g=9.81
PHYSIC=1

pi=np.pi
Pi=np.pi

def sin(x):
	return np.sin(x)
def cos(x):
	return np.cos(x)

def pB_ex(x):
	return 1000.0-100*sin((1/45000.0)*x*Pi)

T0=300
deltaT=50.0

# call functions defined in atm.py
data_f90(dt,g,DEBUG,x0,xf,Nx,p0,Np,PHYSIC)
latm.atm_data.deltat=deltaT
pB=np.zeros((Nx+1),'d',order='f')
for i in range(Nx+1):
	pB[i]=pB_ex(latm.atm_data.x[i,0])
mesh(pB,pA,Nx,Np)

"""*********************************************************
To test mesh, we define our test functions as 
T(x,p) = u(x,p) = w(x,p) = x*p. Later, we can give exact functions
to this functions.
*********************************************************"""

def T_ex(x,p):
	return 300-(1-p/p0)*deltaT

T=np.zeros((Nx+2,Np+2),'d',order='f')
T=T_ex(latm.atm_data.xm,latm.atm_data.pm)

z=fz(latm.atm_data.pm, deltaT, g, p0, T0)

#Sketch fucntion with contour and colorbar.
plt.figure(1)
plt.clf()
plt.contourf(latm.atm_data.xm,z,T)
plt.colorbar(extend="both",format="%.2e")

plt.figure(2)
plt.clf()
plt.contourf(latm.atm_data.xm,latm.atm_data.pm,T)
plt.colorbar(extend="both",format="%.2e")

plt.show()
