#!/usr/bin/env python2.7

# Resolution  with source term of Shallow Water equations with Coriolis force :
# dQ/dt + dF(Q)/dx + dG(Q)/dy + C(Q) = S(t)
#
# With :
# Q = transpose(h,U,V) (data)
#	F = transpose(Fh,Fu,Fv)
#	G = transpose(Gh,Gu,Gv)
#	C = transpose(Ch,Cu,Cv) (coriolis force which equals to f)
#	S = transpose(Sh,SU,SV) (source terme)
#
#	h = height
#	U = h*u
#	V=h*v
#	u = velocity on x-axe
#	v = velocity on y-axe
#
# Method of resolution is fluxes in space and Runge Kutta order 2 in time

import numpy
import matplotlib.pyplot as plt
import libatm as latm
from atm import *
from mlplot import *
import sys
import os
import time


print '######################################'
print '#### Test Phi                     ####'
print '######################################'
#print 'method =', method

dt=1E-3
Nx=400
Np=400
x0=0.0
xf=10.0
p0=10.0
pA=1.0
PHYSIC=0
g=0
DEBUG=0
t0=10.0
nt0=0


# Exact Solution
def sin(x):
	return np.sin(x)
def cos(x):
	return np.cos(x)

pi=np.pi
Pi=np.pi
def pB_ex(x):
	return 10.0

def T_ex(x,p,t):
	return sin(2.0*x*Pi/xf)*(p**2-2.0*pB_ex(x)*p)
def q_ex(x,p,t):
	return sin(2.0*x*Pi/xf)*(p**2-2.0*pB_ex(x)*p)

def u_ex(x,p,t):
	return sin(2.0*x*Pi/xf)*(p**2-2.0*pB_ex(x)*p)

def w_ex(x,p,t):
	return sin(2.0*x*Pi/xf)

def wb(x,t):
	return w_ex(x,pB_ex(x),t)

pB0=p0
R=287.0

	
def phi_dx(x,p,t):
	return -5912.200000-12979.07411* x-7.161972438* sin(0.6283185307 *x)* cos(0.6283185307 *x)+20663.99999 *sin(0.6283185307* x)+5912.200000 *(cos(0.6283185307 *x))**2 - (5596.500000* sin(0.2000000000* x *Pi)+143.5000000 *sin(0.2000000000 *x *Pi) *p**2-5740.0* sin(0.2000000000 *x *Pi)* p)
	
	

# Creation of the mesh
data_f90(dt,g,DEBUG,x0,xf,Nx,p0,Np,PHYSIC)

pB=np.zeros((Nx+1),'d',order='f')
for i in range(Nx+1):
	pB[i]=pB_ex(latm.atm_data.x[i,0])
mesh(pB,pA,Nx,Np)

#print "pm=",latm.atm_data.pm
#print "xm",latm.atm_data.xm
#print "p=",latm.atm_data.p
#print "x",latm.atm_data.x
#print "n",latm.atm_data.normal
#print "ac",latm.atm_data.ac

# Computation
t=t0
latm.atm_data.nt=nt0
X=latm.atm_data.xm
P=latm.atm_data.pm
X2=latm.atm_data.xm2
P2=latm.atm_data.pm2
X3=latm.atm_data.xm3
P3=latm.atm_data.pm3


latm.atm_data.a=np.zeros((2,Nx+2,Np+2),'d',order='f')
latm.atm_data.u=np.zeros((Nx+1,Np+2),'d',order='f')
latm.atm_data.ub=np.zeros((Nx+1),'d',order='f')
latm.atm_data.w=np.zeros((Nx+2,Np+1),'d',order='f')

#print latm.atm_data.dzdx
#print 'dp', latm.atm_data.dp, 'dx', latm.atm_data.dx
#print 'xm',latm.atm_data.xm
#print 'x',latm.atm_data.x
#print 'pm',latm.atm_data.pm
#print 'p',latm.atm_data.p
#print np.size(latm.atm_data.z), np.size(latm.atm_data.xm)

Texac=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
Tbar=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
qexac=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
uexac=np.zeros((Nx+1,Np+2),'d',order='FORTRAN')
wexac=np.zeros((Nx+2,Np+1),'d',order='FORTRAN')

latm.atm_data.sw=np.zeros((Nx+2,Np+1),'d',order='f')
latm.atm_data.s=np.zeros((2,Nx+2,Np+2),'d',order='FORTRAN')
latm.atm_data.sdt=np.zeros((2,Nx+2,Np+2),'d',order='FORTRAN')
latm.atm_data.sdt2=np.zeros((2,Nx+2,Np+2),'d',order='FORTRAN')
latm.atm_data.s2=np.zeros((Nx+1,Np+2),'d',order='FORTRAN')
latm.atm_data.s2dt=np.zeros((Nx+1,Np+2),'d',order='FORTRAN')
latm.atm_data.s2dt2=np.zeros((Nx+1,Np+2),'d',order='FORTRAN')

latm.atm_data.a[0]=T_ex(X,P,t)
latm.atm_data.a[1]=q_ex(X,P,t)
latm.atm_data.u=u_ex(X2,P2,t)
latm.atm_data.w=w_ex(X3,P3,t)
	
uexac=u_ex(X2,P2,t)
latm.atm_data.ub=u_ex(X2[:,Np+1],P2[:,Np+1],t)

#print "w",latm.atm_data.w

print "#### Computation..."



latm.fluxes.borders1()
latm.fluxes.borders2()
latm.fluxes.borders3()
latm.fluxes.res2()
latm.fluxes.fphi()



phi_ex = phi_dx(X,P,t)
errphi=((latm.atm_data.phi[1:Nx+1,1:Np+1]-phi_ex[1:Nx+1,1:Np+1])**2)*latm.atm_data.ac
errphi=np.sqrt(np.sum(errphi))/np.sqrt(np.sum(phi_ex[1:Nx+1,1:Np+1]**2.0*latm.atm_data.ac))
print "L**2 error is " , errphi


print '######################################'
print '#### End Test                     ####'
print '######################################'
