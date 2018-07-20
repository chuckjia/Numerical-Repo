#!/usr/bin/env python2.7

import matplotlib.pyplot as plt
import numpy as np
from atm import *
import libatm as latm

Nx=4
Np=4
DEBUG=0
x0=0.0
xf=10.0
p0=10.0
pA=1.0
PHYSIC=0
dt=1.0
g=0.0
def u(x,p):
	return x
def pB_ex(x):
	return 10.0

data_f90(dt,g,DEBUG,x0,xf,Nx,p0,Np,PHYSIC)
pB=np.zeros((Nx+1),'d',order='f')
for i in range(Nx+1):
	pB[i]=pB_ex(latm.atm_data.x[i,0])
mesh(x0,xf,pB,pA,Nx,Np)
latm.atm_data.bc=1
print "diag=",-2.0/latm.atm_data.dx**2,"low=",1.0/latm.atm_data.dx**2,"up=",1/latm.atm_data.dx**2

LU=latm.fluxes.lu(Nx,1.0)
print "LU(1)",LU[0]
print "LU(2)",LU[1]
print "LU(3)",LU[2]
print "LU(4)",LU[3]
print "LU(5)",LU[4]
LU=latm.fluxes.mlu(Nx,1.0)
print "L", LU[0]
print "U", LU[1]

A=np.zeros((Nx,Nx),'d')
for i in range(Nx):
		for j in range(Nx):
			A[i,j]=np.sum(LU[0,i,:]*LU[1,:,j])
			
print "A", A

latm.atm_data.a=np.zeros((3,Nx+2,Np+2),'d',order='f')
latm.atm_data.lda=np.zeros((Nx),'d',order='f')

latm.atm_data.a[2]=u(latm.atm_data.xm,latm.atm_data.pm)
print "u",latm.atm_data.a[2]	

latm.fluxes.projection(1.0)
print "lambda=",latm.atm_data.lda
	

