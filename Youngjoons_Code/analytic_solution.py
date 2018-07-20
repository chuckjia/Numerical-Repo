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

# Data
SAVEFIG = 1
FIG_NAME = 'analytic_solution'
PLOTFIG = 0 # 1: plot, 0 : no plot
titleT = "T"
titleq = "q"
colors1=np.linspace(-10.5,0,20)
colors2=np.linspace(-0.25,1.75,20)
colors3=[]
colors4=[]
PHYSIC = 1
DEBUG=0
bc=1

# Data for the Mesh
p0=1000.0
#pA=10.0
pA=100.0
x0=0.0
xf=50000.0
Nx=1000#1000#1000#2500
Np=100#100#100#45
u0=0.02
g=9.8
deltaT=50.0
T0=300.0

# Time discretization
nt0=0
t0=0.0
dt=1E-1
tf=t0+1E5*dt
pi=np.pi

def sin(x):
	return np.sin(x)
def cos(x):
	return np.cos(x)
def exp(x):
	return np.exp(x)	

pi=np.pi
Pi=np.pi

def pB_ex(x):
	pbt=np.array(x)
	if len(np.shape(x))==2:
		for i in range(len(x[:,0])):
			for j in range(len(x[0,:])):
				#if (x[i,j]>15000) and (x[i,j]<= 25000):
					##pbt[i,j]=-x[i,j]/10.0+3300
					#pbt[i,j]= 1000 -150*exp(-((x[i,j]-25000)**2)/4500**2)
				#elif (x[i,j]>25000) and (x[i,j]<= 35000):
					##pbt[i,j]=x[i,j]/10.0-2000
					#pbt[i,j]= 1000 -150*exp(-((x[i,j]-25000)**2)/4500**2)
				#else:
					 #pbt[i,j]=1000.0
				pbt[i,j] = 1000 -200*exp(-((x[i,j]-25000)**2)/3000**2)				
	else:
		for i in range(len(x)):
			#if (x[i]>15000) and (x[i]<= 25000):
				##pbt[i]=-x[i]/10.0+3000
				#pbt[i]= 1000 -150*exp(-((x[i]-25000)**2)/4500**2)
			#elif (x[i]>25000) and (x[i]<= 35000):
				##pbt[i]=x[i]/10.0-2000
				#pbt[i]= 1000 -150*exp(-((x[i]-25000)**2)/4500**2)
			#else:
				 #pbt[i]=1000.0
			pbt[i] = 1000.0 -200.0*exp(-((x[i]-25000)**2)/3000**2)						 
	return pbt

# Exact Solution
def T_bar(x,p):
	return T0-(1-p/p0)*deltaT


R=287.0
Cp=1004.0
U0 = 10.0
N = np.sqrt(T0*(1-R/Cp)/p0)
ld = N/U0

def w(x,p):
	return -U0*(-(1/22500)*x+10/9)*exp(-(1/9000000)*(x-25000)**2)*(exp(ld*p)-exp(-ld*p))/(exp(ld*p0)-exp(-ld*p0))

def b(x,p):
	return -.3388528546e-18*exp(-69.44444444+.4628639410e-1*p)+.3388528546e-18*exp(-69.44444444-.4628639410e-1*p)+.3388528546e-18*exp(-.1111111111e-6*x**2+.5555555556e-2*x-69.44444444+.4628639410e-1*p)-.3388528546e-18*exp(-.1111111111e-6*x**2+.5555555556e-2*x-69.44444444-.4628639410e-1*p)
	 
# Creation of the mesh
data_f90(dt,g,DEBUG,x0,xf,Nx,p0,Np,PHYSIC)
latm.atm_data.deltat=deltaT
pB=np.zeros((Nx+1),'d',order='f')
#for i in range(Nx+1):
pB=pB_ex(latm.atm_data.x[:,0])
mesh(x0,xf,pB,pA,Nx,Np)
CellK(Nx,Np,latm.atm_data.x,latm.atm_data.p,latm.atm_data.xm,latm.atm_data.pm)
latm.atm_data.bc=bc


# Computation
t=t0
latm.atm_data.nt=nt0
X=latm.atm_data.xm
P=latm.atm_data.pm
X2=latm.atm_data.xm2
P2=latm.atm_data.pm2
X3=latm.atm_data.xm3
P3=latm.atm_data.pm3

z=fz(P,deltaT, g, p0, T0)
z2=fz(P2,deltaT, g, p0, T0)
z3=fz(P3,deltaT, g, p0, T0)

X=X[:,Np+2-20:]
P=P[:,Np+2-20:]
z=z[:,Np+2-20:]



plot_contour(1,w(X,P),b(X,P),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Analytic $\omega$","Analytic b",FIG_NAME,SAVEFIG,PLOTFIG)
