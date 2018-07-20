#!/usr/bin/env python2.7

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import libatm as latm
from atm import *
from mlplot import *
import sys
import os
import time



Nx=500
Np=500
DEBUG=0
x0=0.0
xf=50000.0
p0=1000.0
pA=100.0
PHYSIC=0
dt=1.0
g=0.0
bc=1

def pB_ex(x):
	pbt=np.array(x)
	if len(np.shape(x))==2:
		for i in range(len(x[:,0])):
			for j in range(len(x[0,:])):
				pbt[i,j] = 1000 -200*exp(-((x[i,j]-25000)**2)/3000**2)				
	else:
		for i in range(len(x)):
			pbt[i] = 1000 -200*exp(-((x[i]-25000)**2)/3000**2)						 
	return pbt
	
def sin(x):
	return np.sin(x)
def cos(x):
	return np.cos(x)
	
def exp(x):
	return np.exp(x)	
	
def u(x,p):
	return cos(2*pi/xf*x)*(p-100)
	#return cos(6*pi/xf*x)*sin(2*pi*p/100)	

pi=np.pi
Pi=np.pi	
	
bc = 1	



data_f90(dt,g,0,x0,xf,Nx,p0,Np,1)

pB=np.zeros((Nx+1),'d',order='f')
#for i in range(Nx+1):
pB=pB_ex(latm.atm_data.x[:,0])
mesh(x0,xf,pB,pA,Nx,Np)
latm.atm_data.bc=bc
dp = latm.atm_data.dp
dx = latm.atm_data.dx
deltaT = 50.0
# Computation
#t=t0
#latm.atm_data.nt=nt0
X=latm.atm_data.xm
P=latm.atm_data.pm
X2=latm.atm_data.xm2
P2=latm.atm_data.pm2
X3=latm.atm_data.xm3
P3=latm.atm_data.pm3
latm.atm_data.zm=np.zeros((Nx+2,Np+2),'d',order='f')
#z=fz(P,deltaT, g, p0, T0)
#z2=fz(P2,deltaT, g, p0, T0)
#z3=fz(P3,deltaT, g, p0, T0)

latm.atm_data.a=np.zeros((3,Nx+2,Np+2),'d',order='f')


latm.atm_data.a[0]=0.0
latm.atm_data.a[1]=0.0
latm.atm_data.a[2]=u(X,P)

#latm.fluxes.Borders2()


#before=CompCond(Nx,Np,dp,dx,latm.atm_data.a[2])

latm.fluxes.projection2(1.0)


#latm.fluxes.Borders2()

#for j in range(Np+2):
#	latm.atm_data.a[2,0 ,j]=latm.atm_data.a[2,1   ,j]
#	latm.atm_data.a[2,Nx+1,j]=latm.atm_data.a[2,Nx,j]
		

after=CompCond(Nx,Np,dp,dx,P,latm.atm_data.a[2])


plot_2d(1,latm.atm_data.ld,after,0,0,latm.atm_data.xm[1:Nx+1,0],latm.atm_data.xm[1:Nx+1,0],"before projection","after projection","test projection1",1,0)

plot_2d(2,latm.atm_data.ld-after,latm.atm_data.ldax,0,0,latm.atm_data.xm[1:Nx+1,0],latm.atm_data.xm[1:Nx+1,0],"Difference before-after","lambda_x","test projection2",1,0)

plot_contour(3,latm.atm_data.a[2],(latm.atm_data.a[2,2:Nx+2,:]-latm.atm_data.a[2,:Nx,:])/(2*latm.atm_data.dx),0,0,X,P,X[1:Nx+1,:],P[1:Nx+1,:],[],[],"u","$u_x$","du",1,0)

plot_contour(4,u(X,P),latm.atm_data.a[2],0,0,X,P,X,P,[],[],"u before","u after","u",1,0)
