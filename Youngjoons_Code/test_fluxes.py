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

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import libatm as latm
from atm import *
from mlplot import *
import sys
import os
import time

os.chdir(sys.argv[1])
sys.path.append(sys.argv[1])

from param import *

DATA_DIR = 'DATA/'
FIG_FILES= "IMG/"+FIG_NAME

if (CLEAN=='yes'):
	if os.path.exists(sys.argv[1]+'/DEBUG')==True:
		for root, dirs, files in os.walk(sys.argv[1]+'/DEBUG', topdown=False):
			for clean_name in files:
				#print root+"/"+clean_name
				os.remove(root+"/"+clean_name)
			for clean_name in dirs:
				os.rmdir(root+"/"+clean_name)
	if os.path.exists(sys.argv[1]+'/DATA')==True:
		for root, dirs, files in os.walk(sys.argv[1]+'/DATA', topdown=False):
			for clean_name in files:
				#print root+"/"+clean_name
				os.remove(root+"/"+clean_name)
			for clean_name in dirs:
				os.rmdir(root+"/"+clean_name)
	if os.path.exists(sys.argv[1]+'/IMG')==True:
		for root, dirs, files in os.walk(sys.argv[1]+'/IMG', topdown=False):
			for clean_name in files:
				#print root+"/"+clean_name
				os.remove(root+"/"+clean_name)
			for clean_name in dirs:
				os.rmdir(root+"/"+clean_name)

# make directory on designate path
if (DEBUG==1 and os.path.exists(sys.argv[1]+'/DEBUG')==False):
		os.mkdir(sys.argv[1]+'/DEBUG')
if (SAVEFIG==1 and os.path.exists(sys.argv[1]+'/IMG')==False):
		os.mkdir(sys.argv[1]+'/IMG')
if (SAVE2TXT_DATA==1 and os.path.exists(sys.argv[1]+'/DATA')==False):
		os.mkdir(sys.argv[1]+'/DATA')

print '######################################'
print '#### Resolution                   ####'
print '######################################'
#print 'method =', method

# Creation of the mesh
data_f90(dt,g,DEBUG,x0,xf,Nx,p0,Np,PHYSIC)
latm.atm_data.deltat=deltaT
pB=np.zeros((Nx+1),'d',order='f')
#for i in range(Nx+1):
pB=pB_ex(latm.atm_data.x[:,0])

mesh(x0,xf,pB,pA,Nx,Np)
CellK(Nx,Np,latm.atm_data.x,latm.atm_data.p,latm.atm_data.xm,latm.atm_data.pm)
CellK2(Nx,Np,latm.atm_data.xm,latm.atm_data.pm)
CellK3(Nx,Np,latm.atm_data.x,latm.atm_data.p,latm.atm_data.xm,latm.atm_data.pm,bc)
#CellK3(Nx,Np,latm.atm_data.xm2,latm.atm_data.pm3,latm.atm_data.xm,latm.atm_data.pm,bc)
coeff_r(Nx,Np,latm.atm_data.xm,latm.atm_data.xm3,latm.atm_data.pm,latm.atm_data.pm2)
latm.atm_data.bc=bc
latm.atm_data.th1=th1
latm.atm_data.th2=th2
latm.atm_data.scheme=scheme

# Computation
t=t0
latm.atm_data.nt=nt0
X=latm.atm_data.xm
P=latm.atm_data.pm
X2=latm.atm_data.xm2
P2=latm.atm_data.pm2
X3=latm.atm_data.xm3
P3=latm.atm_data.pm3


#latm.atm_data.coeff3=0.25*np.ones(np.shape(latm.atm_data.coeff3))


z=fz(P,deltaT, g, p0, T0)
latm.atm_data.zm=np.zeros((Nx+2,Np+2),'d',order='f')
latm.atm_data.zm=z
z2=fz(P2,deltaT, g, p0, T0)
z3=fz(P3,deltaT, g, p0, T0)

latm.atm_data.a=np.zeros((3,Nx+2,Np+2),'d',order='f')
latm.atm_data.w=np.zeros((Nx+2,Np+2),'d',order='f')


Texac=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
Tbar=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
qexac=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
uexac=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
wexac=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
temp =np.zeros((Nx+2,Np+2),'d',order='FORTRAN')

latm.atm_data.s=np.zeros((3,Nx+2,Np+2),'d',order='FORTRAN')
latm.atm_data.sdt=np.zeros((3,Nx+2,Np+2),'d',order='FORTRAN')
latm.atm_data.sdt2=np.zeros((3,Nx+2,Np+2),'d',order='FORTRAN')

if (PLOT_ERROR==1):
	ind_t=0
	print tf, PERIODE_ERROR, tf/dt/PERIODE_ERROR
	errT_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
	errq_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
	errU_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
	errw_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
	t_t=np.zeros((tf/dt/PERIODE_ERROR),'d')

if (SOURCE==0):
	latm.atm_data.a[0]=T_init(X,P)
	Tbar=T_bar(X,P)
	latm.atm_data.a[1]=q_init(X,P)
	latm.atm_data.a[2]=u_init(X,P)
	latm.atm_data.w=w_init(X,P)
else :
	latm.atm_data.a[0]=T_ex(X,P,t)
	Texac=T_ex(X,P,t)
	latm.atm_data.a[1]=q_ex(X,P,t)
	qexac=q_ex(X,P,t)
	latm.atm_data.sdt[0]=ST(X,P,t)
	latm.atm_data.sdt[1]=Sq(X,P,t)
	latm.atm_data.a[2]=u_ex(X,P,t)
	uexac=u_ex(X,P,t)
	latm.atm_data.sdt[2]=Su(X,P,t)
	latm.atm_data.w=w_ex(X,P,t)
	wexac=w_ex(X,P,t)



print "#### Computation..."

latm.fluxes.borders1()
	
if (scheme==1):
	latm.central_upwind.dfdp()
	latm.central_upwind.dgdx()
else:
	latm.central_upwind.dfdp2()
	latm.central_upwind.dgdx2()


latm.fluxes.fomega()
latm.fluxes.borders3()



wexac=w_ex(X,P,t)
ac=latm.atm_data.ac
w=latm.atm_data.w

errw=((w[1:Nx+1,1:Np+1]-wexac[1:Nx+1,1:Np+1])**2)*ac
errw=np.sqrt(np.sum(errw))/np.sqrt(np.sum(wexac[1:Nx+1,1:Np+1]**2.0*ac))
print  "log(errw)=",np.log10(errw)


Texac=STx(X,P,t)
qexac=Sqx(X,P,t)
uexac=Sux(X,P,t)

T = latm.atm_data.fluxhori[0]
q = latm.atm_data.fluxhori[1]
u = latm.atm_data.fluxhori[2]

errT=((T-Texac[1:Nx+1,1:Np+1])**2)*ac
errq=((q-qexac[1:Nx+1,1:Np+1])**2)*ac
erru=((u-uexac[1:Nx+1,1:Np+1])**2)*ac
errT=np.sqrt(np.sum(errT))/np.sqrt(np.sum(Texac[1:Nx+1,1:Np+1]**2.0*ac))
errq=np.sqrt(np.sum(errq))/np.sqrt(np.sum(qexac[1:Nx+1,1:Np+1]**2.0*ac))
erru=np.sqrt(np.sum(erru))/np.sqrt(np.sum(uexac[1:Nx+1,1:Np+1]**2.0*ac))

print " log(errT)=",np.log10(errT)," log(errq)=",np.log10(errq)," log(erru)=",np.log10(erru)


Texac=STp(X,P,t)
qexac=Sqp(X,P,t)
uexac=Sup(X,P,t)
T = latm.atm_data.fluxvert[0]
q = latm.atm_data.fluxvert[1]
u = latm.atm_data.fluxvert[2]


errT=((T-Texac[1:Nx+1,1:Np+1])**2)*ac
errq=((q-qexac[1:Nx+1,1:Np+1])**2)*ac
erru=((u-uexac[1:Nx+1,1:Np+1])**2)*ac
errT=np.sqrt(np.sum(errT))/np.sqrt(np.sum(Texac[1:Nx+1,1:Np+1]**2.0*ac))
errq=np.sqrt(np.sum(errq))/np.sqrt(np.sum(qexac[1:Nx+1,1:Np+1]**2.0*ac))
erru=np.sqrt(np.sum(erru))/np.sqrt(np.sum(uexac[1:Nx+1,1:Np+1]**2.0*ac))

print " log(errT)=",np.log10(errT)," log(errq)=",np.log10(errq)," log(erru)=",np.log10(erru)

Texac=ST(X,P,t)
qexac=Sq(X,P,t)
uexac=Su(X,P,t)
T = latm.atm_data.fluxvert[0]+latm.atm_data.fluxhori[0]
q = latm.atm_data.fluxvert[1]+latm.atm_data.fluxhori[1]
u = latm.atm_data.fluxvert[2]+latm.atm_data.fluxhori[2]

errT=((T-Texac[1:Nx+1,1:Np+1])**2)*ac
errq=((q-qexac[1:Nx+1,1:Np+1])**2)*ac
erru=((u-uexac[1:Nx+1,1:Np+1])**2)*ac
errT=np.sqrt(np.sum(errT))/np.sqrt(np.sum(Texac[1:Nx+1,1:Np+1]**2.0*ac))
errq=np.sqrt(np.sum(erru))/np.sqrt(np.sum(qexac[1:Nx+1,1:Np+1]**2.0*ac))
erru=np.sqrt(np.sum(erru))/np.sqrt(np.sum(uexac[1:Nx+1,1:Np+1]**2.0*ac))

print " log(errT)=",np.log10(errT)," log(errq)=",np.log10(errq)," log(erru)=",np.log10(erru)


print '######################################'
print '#### End Resolution               ####'
print '######################################'
