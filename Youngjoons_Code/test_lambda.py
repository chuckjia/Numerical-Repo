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
print '#### Resolution                   ####'
print '######################################'
#print 'method =', method

dt=1E-3
Nx=50
Np=50
x0=0.0
xf=10.0
p0=10.0
pA=1.0
PHYSIC=0
g=0
DEBUG=0
t0=10.0
nt0=0
dt=1E-4


# Exact Solution
def sin(x):
	return np.sin(x)
def cos(x):
	return np.cos(x)

pi=np.pi
Pi=np.pi

def T_ex(x,p,t):
	return cos(x*pi/xf)
def q_ex(x,p,t):
	return x

def u_ex(x,p,t):
	return cos(t*pi)*sin(x*2.0*pi/xf)*(p**2-2.0*pB_ex(x)*p)

def wb(x,t):
	return w_ex(x,pB_ex(x),t)

pB0=p0
R=287.0

def pB_ex(x):
	return 10.0

def ld_exdxx(x,t):
	return 1/(pB0-pA)*cos(t*pi)*cos(x*2.0*pi/xf)*(2.0*pi/xf)*(-2.0*pB0**3/3.0 - (pA**3/3.0-pB0*pA**2))
	
def ld_exdx(x,t):
	return 1/(pB0-pA)*cos(t*pi)*sin(x*2.0*pi/xf)*(-2.0*pB0**3/3.0 - (pA**3/3.0-pB0*pA**2))
	
def ld_ex1(x,t):
	return -1/(pB0-pA)*cos(t*pi)*cos(x*2.0*pi/xf)*xf/(2.0*pi)*(-2.0*pB0**3/3.0 - (pA**3/3.0-pB0*pA**2))
	
def ld_ex(x,t):
	return -1/(pB0-pA)*cos(t*pi)*cos(x*2.0*pi/xf)*xf/(2.0*pi)*(-2.0*pB0**3/3.0 - (pA**3/3.0-pB0*pA**2)) - ld_ex1(xf-(xf-x0)/float(Nx),t)
	
def pu_ex(x,p,t):
	return -1/(pB0-pA)*cos(t*pi)*sin(x*2.0*pi/xf)*(-2.0*pB0**3/3.0 - (pA**3/3.0-pB0*pA**2))+cos(t*pi)*sin(x*pi/xf)*(p**2-2.0**pB0*p);
	
def phi_dx(x,p,t):
	return (.5000000000*(104.9760000*cos(t*Pi)**2*sin(.1000000000*x*Pi)**2*Pi**2\
		+2.824605956*10**(-8)*cos(3.141592654*t)*cos(.3141592654*x)*cos(t*Pi)*\
		cos(.1000000000*x*Pi)*Pi-8.873761320*10**(-8)*cos(3.141592654*t)*sin(.3141592654*x)\
		*cos(t*Pi)*sin(.1000000000*x*Pi)-104.9760000*cos(t*Pi)**2*cos(.1000000000*x*Pi)**2*Pi**2\
		+.4050000000*R*cos(.1000000000*x*Pi)*Pi**2))*x**2+x**2*(-104.9760000*cos(t*Pi)**2*sin(.1000000000*x*Pi)**2\
		*Pi**2-2.824605956*10**(-8)*cos(3.141592654*t)*cos(.3141592654*x)*cos(t*Pi)*cos(.1000000000*x*Pi)*\
		Pi+8.873761320*10**(-8)*cos(3.141592654*t)*sin(.3141592654*x)*cos(t*Pi)*sin(.1000000000*x*Pi)+\
		104.9760000*cos(t*Pi)**2*cos(.1000000000*x*Pi)**2*Pi**2-.4050000000*R*cos(.1000000000*x*Pi)*Pi**2)
	
def w_ex(x,p,t):
	return -.3333333333*cos(t*Pi)*cos(x*Pi/xf)*Pi*(p**3-1.*pA**3)/xf+1.*cos(t*Pi)*\
		cos(x*Pi/xf)*Pi*pB0*(p**2-1.*pA**2)/xf+1.*cos(t*Pi)*cos(x*Pi/xf)*Pi*\
		(-.6666666667*pB0**3-.3333333333*pA**3+pB0*pA**2)*(p-1.*pA)/((pB0-1.*pA)*xf)
	

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
latm.atm_data.ub]u_ex(X2[:,0],P2[:,0],t-dt)
	
uexac=pu_ex(X2,P2,t)

#print "w",latm.atm_data.w

print "#### Computation..."

#print "dp, dx", latm.atm_data.dp, latm.atm_data.dx
#print "hhori", latm.atm_data.hhori
#print "x, xm", latm.atm_data.x, latm.atm_data.xm
#print "p, pm", latm.atm_data.p, latm.atm_data.pm
#print "Cp", latm.atm_data.cp
#print "w(0)", latm.atm_data.w[:,0]
#print "q=",latm.atm_data.a[1]

print latm.atm_data.xm[Nx-1,1]



#latm.fluxes.borders2()
latm.fluxes.projection2()
erru=((latm.atm_data.u[:Nx,1:Np+1]-uexac[:Nx,1:Np+1])**2)*latm.atm_data.ac2
erru=np.sqrt(np.sum(erru))/np.sqrt(np.sum(uexac[:Nx,1:Np+1]**2.0*latm.atm_data.ac2))

print erru

print len(latm.atm_data.ldax),len(X2[:Nx,1]),len(ld_exdx(X2[:Nx,1],t))

print X[1:Nx+1,1]
print "lambda_dxdx",ld_exdxx(X[1:Nx+1,1],t)
print X[0,1], X[Nx+1,1], ld_exdxx(X[0,1],t), ld_exdxx(X[Nx+1,1],t)

print latm.atm_data.ldax/ld_exdx(X2[:Nx,1],t)


plt.figure(1)
plt.clf()
plt.plot(X2[:Nx,1],-latm.atm_data.lda+ld_ex(X[1:Nx+1,1],t))

plt.figure(2)
plt.clf()
plt.plot(X2[:Nx,1],latm.atm_data.lda, 'r')
plt.plot(X2[:Nx,1],ld_ex(X[1:Nx+1,1],t),'g')

plt.figure(3)
plt.clf()
plt.plot(X2[:Nx,1],latm.atm_data.ldax-ld_exdx(X2[:Nx,1],t))


plt.show()

print '######################################'
print '#### End Resolution               ####'
print '######################################'
