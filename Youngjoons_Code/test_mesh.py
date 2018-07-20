#!/usr/bin/env python2.7

import numpy
import libatm as latm
from atm import *
from mlplot import *
import sys
import os
import time
import matplotlib.pyplot as plt

Nx=10
Np=10
DEBUG=0
x0=0.0
xf=10.0
p0=10.0
pA=0.0
PHYSIC=0
dt=1.0
g=0.0
bc=1

def pB_ex(x):
	pbt=np.array(x)
	if len(np.shape(x))==2:
		for i in range(len(x[:,0])):
			for j in range(len(x[0,:])):
				pbt[i,j] = 10.0 -4.0*exp(-((x[i,j]-5.0)**2))				
				
	else:
		for i in range(len(x)):
			pbt[i] = 10.0 -4.0*exp(-((x[i]-5.0)**2))				
	return pbt
	
def sin(x):
	return np.sin(x)
def cos(x):
	return np.cos(x)
	
def exp(x):
	return np.exp(x)	
	
def u(x,p):
	return cos(2*pi/xf*x)*(p**2+p)
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
X=latm.atm_data.x
P=latm.atm_data.p



fig = plt.figure()
ax = fig.add_subplot(111)

coord = np.zeros((Nx+1,2,Np+1))
temp_coord = np.zeros((2,Np+1))

for i in range(Nx+1):
	for j in range(Np+1):
		coord[i,0,j] = X[i,j]
		coord[i,1,j] = P[i,j]

for i in range(Nx+1):
	temp_coord = coord[i,:,:]
	temp = zip(*temp_coord)
	xs,ys = zip(*temp)
	ax.plot(xs, ys,"k")

for i in range(Nx+1):
	for j in range(Np+1):
		coord[i,0,j] = X[j,i]
		coord[i,1,j] = P[j,i]


for i in range(Nx+1):
	temp_coord = coord[i,:,:]
	temp = zip(*temp_coord)
	xs,ys = zip(*temp)
	ax.plot(xs, ys,"k")


ax.set_xlim(-2, 12)
ax.set_ylim(-2, 12)
ax.set_title('Example of the grids')
ax.set_xlabel('x-axis')
ax.set_ylabel('p-axis')
plt.show()


#CellK3(Nx,Np,X2,P3,X,P,0)
#CellK2(Nx,Np,latm.atm_data.xm,latm.atm_data.pm)
#CellK3(Nx,Np,latm.atm_data.x,latm.atm_data.p,X,P,0)
#print 'X = ',X
#print 'P = ',P
#print 'X2 = ',X2
#print 'P2 = ',P2
#print 'X3 = ',X3
#print 'P3 = ',P3
#print 'coeff3(1) = ',latm.atm_data.coeff3[:,:,0]
#print 'coeff3(2) = ',latm.atm_data.coeff3[:,:,1]
#print 'coeff3(3) = ',latm.atm_data.coeff3[:,:,2]
#print 'coeff3(4) = ',latm.atm_data.coeff3[:,:,3]

#print 'eqaul to 1', latm.atm_data.coeff3[:,:,0]+latm.atm_data.coeff3[:,:,1]+latm.atm_data.coeff3[:,:,2]+latm.atm_data.coeff3[:,:,3]

#print 'm3(1) = ',latm.atm_data.m3[:,:,0]
#print 'm3(2) = ',latm.atm_data.m3[:,:,1]
#print 'm3(3) = ',latm.atm_data.m3[:,:,2]
#print 'm3(4) = ',latm.atm_data.m3[:,:,3]


#print 'm2(1) = ',latm.atm_data.m2[:,:,0]
#print 'm2(2) = ',latm.atm_data.m2[:,:,1]
#print 'm2(3) = ',latm.atm_data.m2[:,:,2]
#print 'm2(4) = ',latm.atm_data.m2[:,:,3]




#latm.atm_data.a=np.zeros((3,Nx+2,Np+2),'d',order='f')


#latm.atm_data.a[0]=0.0
#latm.atm_data.a[1]=0.0
#latm.atm_data.a[2]=u(X,P)

#latm.fluxes.Borders2()


#before=CompCond(Nx,Np,dp,dx,latm.atm_data.a[2])

#latm.fluxes.projection2(1)


#latm.fluxes.Borders2()

#for j in range(Np+2):
#	latm.atm_data.a[2,0 ,j]=latm.atm_data.a[2,1   ,j]
#	latm.atm_data.a[2,Nx+1,j]=latm.atm_data.a[2,Nx,j]
		

#after=CompCond(Nx,Np,dp,dx,latm.atm_data.a[2])


#plot_2d(1,before,after,0,0,latm.atm_data.xm[1:Nx+1,0],latm.atm_data.xm[1:Nx+1,0],"before proje","after proj","test projection1",1,1)

#plot_2d(3,latm.atm_data.ld,latm.atm_data.ldax,0,0,latm.atm_data.xm[1:Nx+1,0],latm.atm_data.xm[1:Nx+1,0],"ld","lambda_x","test projection2",1,1)

