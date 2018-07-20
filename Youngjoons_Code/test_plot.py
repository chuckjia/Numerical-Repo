#!/usr/bin/env python2.7

import matplotlib.pyplot as plt
import numpy as np
from atm import *
import libatm as latm


"""**************************************************************
To call the function data_f90 and mesh, we just define following 
tempoaray data. 
**************************************************************"""
Nx=3
Np=3
#pB=[10.0,8.0,7.0,8.0,10.0]
pB=[10.0,10.0,10.0,10.0]
x0=0.0
xf=10.0
pA=0.0
p0=10.0

theta=0
bc=[0,0]
bc2=[0,0]
DEBUG=0
dt=1
g=1
PHYSIC=0


# call functions defined in atm.py
data_f90(dt,g,DEBUG,x0,xf,Nx,p0,Np,PHYSIC)
mesh(x0,xf,pB,pA,Nx,Np)
x=latm.atm_data.x
p=latm.atm_data.p
xm=latm.atm_data.xm
pm=latm.atm_data.pm

print x
print xm

CellK(Nx,Np,x,p,xm,pm)
print 'm',latm.atm_data.m
print 'coeff',latm.atm_data.coeff


"""*********************************************************
To test mesh, we define our test functions as 
T(x,p) = u(x,p) = w(x,p) = x*p. Later, we can give exact functions
to this functions.
*********************************************************"""

def T_ex(x,p):
	return x*p

def u_ex(x,p):
	return x*p

def w_ex(x,p):
	return x*p


T=np.zeros((Nx+2,Np+2),'d',order='f')
T=T_ex(latm.atm_data.xm,latm.atm_data.pm)
		
u=np.zeros((Nx+1,Np+2),'d',order='f')
u=u_ex(latm.atm_data.xm2,latm.atm_data.pm2)

w=np.zeros((Nx+2,Np+1),'d',order='f')
w=u_ex(latm.atm_data.xm3,latm.atm_data.pm3)


# Print out and check our data for meshes.
#print "dp=",latm.atm_data.dp
#print "x=",latm.atm_data.x
print "xm=",latm.atm_data.xm
#print "p=",latm.atm_data.p
#print "xm_new" ,latm.atm_data.xm_new
print "pm=",latm.atm_data.pm
#print "pm_new" ,latm.atm_data.pm_new
#print "diff_xm",(latm.atm_data.xm_new - latm.atm_data.xm)
#print "diff_pm",(latm.atm_data.pm_new - latm.atm_data.pm)

#print "nx=",latm.atm_data.normal[0]
#print "ny=",latm.atm_data.normal[1]
#print "h horizontal=", latm.atm_data.hhori
#print "ac=", latm.atm_data.ac
#print "ac2=", latm.atm_data.ac2
#print "ac3=", latm.atm_data.ac3

#print "xm2=",latm.atm_data.xm2
#print "pm2=",latm.atm_data.pm2
#print "xm3=",latm.atm_data.xm3
#print "pm3=",latm.atm_data.pm3


"""


np.reshape(latm.atm_data.x[1:,:]-latm.atm_data.x[:Nx,:],Nx*(Np+1))
print "n * (x,y) =", latm.atm_data.normal[0]*np.reshape(latm.atm_data.x[1:,:]-latm.atm_data.x[:Nx,:],Nx*(Np+1),order='F') \
	+latm.atm_data.normal[1]*np.reshape(latm.atm_data.p[1:,:]-latm.atm_data.p[:Nx,:],Nx*(Np+1),order='F')

print "n ^2=", latm.atm_data.normal[0]**2+latm.atm_data.normal[1]**2


#Sketch fucntion with contour and colorbar.
#plt.figure(fig,figsize=(10,5))
plt.figure(1)
plt.clf()
plt.contourf(latm.atm_data.xm,latm.atm_data.pm,T)
plt.colorbar(extend="both",format="%.2e")
plt.ylim(0,10.0)
#plt.show()

#plt.figure(fig,figsize=(10,5))
plt.figure(2)
plt.clf()
plt.contourf(latm.atm_data.xm2,latm.atm_data.pm2,u)
plt.colorbar(extend="both",format="%.2e")
plt.ylim(0,10.0)
#plt.show()
#plt.figure(fig,figsize=(10,5))

plt.figure(3)
plt.clf()
plt.contourf(latm.atm_data.xm3,latm.atm_data.pm3,w)
plt.colorbar(extend="both",format="%.2e")
plt.ylim(0,10.0)

plt.show()

"""
