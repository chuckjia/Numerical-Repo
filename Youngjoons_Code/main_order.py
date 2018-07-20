#!/usr/bin/env python2.6

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

os.chdir(sys.argv[1])
sys.path.append(sys.argv[1])

from param import *

DATA_DIR = 'DATA/'
FIG_FILES= "IMG/"+FIG_NAME

if (CLEAN=='yes'):
	if os.path.exists(sys.argv[1]+'/DEBUG')==True:
		for root, dirs, files in os.walk(sys.argv[1]+'/DEBUG', topdown=False):
			for clean_name in files:
				print root+"/"+clean_name
				os.remove(root+"/"+clean_name)
			for clean_name in dirs:
				os.rmdir(root+"/"+clean_name)
	if os.path.exists(sys.argv[1]+'/DATA')==True:
		for root, dirs, files in os.walk(sys.argv[1]+'/DATA', topdown=False):
			for clean_name in files:
				print root+"/"+clean_name
				os.remove(root+"/"+clean_name)
			for clean_name in dirs:
				os.rmdir(root+"/"+clean_name)
	if os.path.exists(sys.argv[1]+'/IMG')==True:
		for root, dirs, files in os.walk(sys.argv[1]+'/IMG', topdown=False):
			for clean_name in files:
				print root+"/"+clean_name
				os.remove(root+"/"+clean_name)
			for clean_name in dirs:
				os.rmdir(root+"/"+clean_name)

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

it=len(NxT)
logdp=np.zeros(it)
logdx=np.zeros(it)
logerrT=np.zeros(it)
logerrq=np.zeros(it)
logerru=np.zeros(it)

for c in range(it):
	Nx=NxT[c]
	Np=NpT[c]
	
	# Creation of the mesh
	[X,P]=data_f90(dt,theta,u0,g,DEBUG,x0,xf,Nx,p0,pB,Np)

	# Computation
	t=t0
	latm.atm_data.nt=nt0

	latm.atm_data.a=np.zeros((2,Nx,Np),'d',order='f')
	latm.atm_data.u=np.zeros((Nx+1,Np),'d',order='f')
	latm.atm_data.w=np.zeros((Nx,Np+1),'d',order='f')
	latm.atm_data.z=np.zeros((Nx),'d',order='f')
	latm.atm_data.dzdx=np.zeros((Nx),'d',order='f')
	latm.atm_data.dpdx=np.zeros((Nx),'d',order='f')

	Texac=np.zeros((Nx,Np),'d',order='FORTRAN')
	qexac=np.zeros((Nx,Np),'d',order='FORTRAN')
	uexac=np.zeros((Nx+1,Np),'d',order='FORTRAN')
	wexac=np.zeros((Nx,Np+1),'d',order='FORTRAN')

	for i in range(Nx):
		latm.atm_data.z[i]=z(latm.atm_data.xm[i])
		latm.atm_data.dzdx[i]=dzdx(latm.atm_data.xm[i])
		latm.atm_data.dpdx[i]=dpdx(latm.atm_data.xm[i])

	if (PLOT_ERROR==1):
		ind_t=0
		print tf, PERIODE_ERROR, tf/dt/PERIODE_ERROR
		errT_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
		errq_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
		errU_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
		errw_t=np.zeros((tf/dt/PERIODE_ERROR),'d')
		t_t=np.zeros((tf/dt/PERIODE_ERROR),'d')

	if (SOURCE==0):
		for i in range(Nx):
			for j in range(Np):
				latm.atm_data.a[0,i,j]=T_init(latm.atm_data.xm[i],latm.atm_data.pm[j])
				latm.atm_data.a[1,i,j]=q_init(latm.atm_data.xm[i],latm.atm_data.pm[j])
				latm.atm_data.u[i,j]=u_init(latm.atm_data.x[i],latm.atm_data.pm[j])
				latm.atm_data.w[i,j]=w_init(latm.atm_data.xm[i],latm.atm_data.p[j])
		for j in range(Np):
			latm.atm_data.u[Nx+1,j]=u_init(latm.atm_data.x[Nx+1],latm.atm_data.pm[j])
		for i in range(Nx):
			latm.atm_data.w[i,Np+2]=w_init(latm.atm_data.xm[i],latm.atm_data.p[Np+1])

			#latm.atm_data.sdt[0,i]=latm.atm_data.s[0,i]
			#latm.atm_data.sdt2[0,i]=latm.atm_data.s[0,i]

	else :
		latm.atm_data.sw=np.zeros((Nx,Np+1),'d',order='f')
		latm.atm_data.s=np.zeros((2,Nx,Np),'d',order='FORTRAN')
		latm.atm_data.sdt=np.zeros((2,Nx,Np),'d',order='FORTRAN')
		latm.atm_data.sdt2=np.zeros((2,Nx,Np),'d',order='FORTRAN')
		latm.atm_data.s2=np.zeros((Nx,Np),'d',order='FORTRAN')
		latm.atm_data.s2dt=np.zeros((Nx,Np),'d',order='FORTRAN')
		latm.atm_data.s2dt2=np.zeros((Nx,Np),'d',order='FORTRAN')
		for i in range(Nx):
			for j in range(Np):
				latm.atm_data.a[0,i,j]=T_ex(latm.atm_data.xm[i],latm.atm_data.pm[j],t)
				latm.atm_data.a[1,i,j]=q_ex(latm.atm_data.xm[i],latm.atm_data.pm[j],t)
				latm.atm_data.sdt[0,i,j]=ST(latm.atm_data.xm[i],latm.atm_data.pm[j],t)
				latm.atm_data.sdt[1,i,j]=Sq(latm.atm_data.xm[i],latm.atm_data.pm[j],t)
				latm.atm_data.u[i,j]=u_ex(latm.atm_data.x[i],latm.atm_data.pm[j],t)
				latm.atm_data.s2dt[i,j]=Su(latm.atm_data.x[i],latm.atm_data.pm[j],t)
				latm.atm_data.w[i,j]=w_ex(latm.atm_data.xm[i],latm.atm_data.p[j],t)
				latm.atm_data.sw[i,j]=Sw(latm.atm_data.xm[i],latm.atm_data.p[j],t)
		for j in range(Np):
			latm.atm_data.u[Nx,j]=u_ex(latm.atm_data.x[Nx],latm.atm_data.pm[j],t)
		for i in range(Nx):
			latm.atm_data.w[i,Np]=w_ex(latm.atm_data.xm[i],latm.atm_data.p[Np],t)
			latm.atm_data.sw[i,Np]=Sw(latm.atm_data.xm[i],latm.atm_data.p[Np],t)

	print "#### Computation..."

	t_before=time.time()

	while t<tf-dt/2.0:

		#print "t=",t

		if (SOURCE==1):
			for i in range(Nx):
				for j in range(Np):
					latm.atm_data.s[0,i,j]=latm.atm_data.sdt[0,i,j]
					latm.atm_data.s[1,i,j]=latm.atm_data.sdt[1,i,j]
					latm.atm_data.s2[i,j]=latm.atm_data.s2dt[i,j]
					latm.atm_data.sdt[0,i,j]=ST(latm.atm_data.xm[i],latm.atm_data.pm[j],t+dt)
					latm.atm_data.sdt[1,i,j]=Sq(latm.atm_data.xm[i],latm.atm_data.pm[j],t+dt)
					latm.atm_data.s2dt[i,j]=Su(latm.atm_data.x[i],latm.atm_data.pm[j],t+dt)	
					latm.atm_data.sdt2[0,i,j]=ST(latm.atm_data.xm[i],latm.atm_data.pm[j],t+dt/2.0)
					latm.atm_data.sdt2[1,i,j]=Sq(latm.atm_data.xm[i],latm.atm_data.pm[j],t+dt/2.0)
					latm.atm_data.s2dt2[i,j]=Su(latm.atm_data.x[i],latm.atm_data.pm[j],t+dt/2.0)
					latm.atm_data.sw[i,j]=Sw(latm.atm_data.xm[i],latm.atm_data.p[j],t)
			for i in range (Nx):
				latm.atm_data.sw[i,Np]=Sw(latm.atm_data.xm[i],latm.atm_data.p[Np],t)


		if (time_method==2):
			RK2Source()
		if (time_method==4):
			RK4Source()

		#ATTENTION POUR CAS TEST:
		#for i in range(Nx):
			#for j in range(Np+1):
				#latm.atm_data.w[i,j]=w_ex(latm.atm_data.xm[i],latm.atm_data.p[j],t)

		t=t+dt
		latm.atm_data.nt=latm.atm_data.nt+1

		if (SOURCE==1):
			for i in range(Nx):
				for j in range(Np):
					Texac[i,j]=T_ex(latm.atm_data.xm[i],latm.atm_data.pm[j],t)
					qexac[i,j]=q_ex(latm.atm_data.xm[i],latm.atm_data.pm[j],t)
					uexac[i,j]=u_ex(latm.atm_data.x[i],latm.atm_data.pm[j],t)
					wexac[i,j]=w_ex(latm.atm_data.xm[i],latm.atm_data.p[j],t)
			for j in range(Np):
				uexac[Nx,j]=u_ex(latm.atm_data.x[Nx],latm.atm_data.pm[j],t)
			for i in range(Nx):
				wexac[i,Np]=w_ex(latm.atm_data.xm[i],latm.atm_data.p[Np],t)

			[errT,errq,erru,errw]=L2_norm(Nx,Np,latm.atm_data.dx,latm.atm_data.dp,latm.atm_data.a[0],latm.atm_data.a[1],latm.atm_data.u,latm.atm_data.w,Texac,qexac,uexac,wexac)
			print "t=",t," log(dx)=",numpy.log10(latm.atm_data.dx)," log(dp)=",numpy.log10(latm.atm_data.dp)," log(errT)=",numpy.log10(errT)," log(errq)=",numpy.log10(errq)," log(erru)=",numpy.log10(erru)," log(errw)=",numpy.log10(errw)
		else:
			 print "t=",t," norm L2 of T=",np.sqrt(latm.atm_data.dp*np.inner(latm.atm_data.a[0],latm.atm_data.a[0])), " norme L2 of q=",np.sqrt(latm.atm_data.dp*np.inner(latm.atm_data.a[1],latm.atm_data.a[1]))


	t_after=time.time()

	t_elapsed=t_after-t_before

	print "#### Elapsed time = ",t_elapsed

	logdx[c]=numpy.log10(latm.atm_data.dx)
	logdp[c]=numpy.log10(latm.atm_data.dp)
	logerrT[c]=numpy.log10(errT)
	logerrq[c]=numpy.log10(errq)
	logerru[c]=numpy.log10(erru)

	print '######################################'
	print '#### End Resolution               ####'
	print '######################################'

plt.figure(1)	
plt.clf()
plt.title('Error in space a t='+str(t))
plt.plot(logdx,logerrT,'b-',logdx,logerrq,'r-',logdx,logerru,'g-')
filename="IMG/"+"Error_in_space"
plt.legend( ('log(errT)', 'log(errq)', 'log(erru)') )
plt.savefig(filename)

print "Slopes..."
print "T :", (logerrT[it-1]-logerrT[it-2])/(logdx[it-1]-logdx[it-2])
print "q :", (logerrq[it-1]-logerrq[it-2])/(logdx[it-1]-logdx[it-2])
print "u :", (logerru[it-1]-logerru[it-2])/(logdx[it-1]-logdx[it-2])
