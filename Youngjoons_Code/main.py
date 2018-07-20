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
import struct

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
latm.atm_data.ep_dump=Ep_dump

#latm.atm_data.rs=1/2.0*np.ones(np.shape(latm.atm_data.rs))
#latm.atm_data.rw=1/2.0*np.ones(np.shape(latm.atm_data.rs))

# Computation
t=t0
latm.atm_data.nt=nt0
X=latm.atm_data.xm
P=latm.atm_data.pm
X2=latm.atm_data.xm2
P2=latm.atm_data.pm2
X3=latm.atm_data.xm3
P3=latm.atm_data.pm3


#print latm.atm_data.normal
#print np.max(X[:,2:Np+2]-X[:,0:Np]),latm.atm_data.dx


z=fz(P,deltaT, g, p0, T0)
latm.atm_data.zm=np.zeros((Nx+2,Np+2),'d',order='f')
latm.atm_data.zm=z
z2=fz(P2,deltaT, g, p0, T0)
z3=fz(P3,deltaT, g, p0, T0)

latm.atm_data.a=np.zeros((3,Nx+2,Np+2),'d',order='f')
latm.atm_data.w=np.zeros((Nx+2,Np+2),'d',order='f')

#print "xm=",X[Nx/2,:], X[Nx/2-5,:], X[Nx/2+5,:]
#print "pm=",P#,"\n pB_ex=",pB_ex(X)
#print "xm", X
#print "zm=",z
#print_mesh(Nx,Np,P)
#print latm.atm_data.x[:,0], pB

Texac=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
latm.atm_data.tbar=np.zeros((Nx+2,Np+2),'d',order='FORTRAN')
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


# RANDOM q_s
random_qs(ep,Nx,Np) 
Average_Sol(Nx,Np)

#if ((t % 2*dt) == 0):
	#print 'average=',t
	#Average_Sol(Nx,Np)
	
	


if (SOURCE==0):
	latm.atm_data.a[0]=T_init(X,P)
	latm.atm_data.tbar=T_bar(X,P)
	latm.atm_data.a[1]=q_init(X,P)+ latm.atm_data.rqs
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


#for i in range(Nx+2):
	#for j in range(Np+2):
		#if latm.atm_data.a[0,i,j] < 275.0:
			#latm.atm_data.a[0,i,j] = 275.0

temp_qs = q_s(X,P)
# Computation of initial values
if PHYSIC==1:
	latm.fluxes.projection2(1.0)
	#latm.fluxes.fphix()
	#latm.fluxes.fomega()
	
elif PHYSIC==2:
	latm.fluxes.projection2(1.0)
	latm.fluxes.fphix()
	#latm.fluxes.fomega()

latm.fluxes.fomega()

print "t=",t


#print 'phix = ', latm.atm_data.phix

if (SAVE2TXT_DATA==1):
	save_txt(DATA_DIR,t,Nx,Np,X,P,latm.atm_data.a[0],latm.atm_data.a[1],latm.atm_data.a[2],latm.atm_data.w)

latm.atm_data.a[2,0,:]=latm.atm_data.a[2,1,:]
#FILE_ID =  open('test.dat', 'w')    
#FILE_ID.write(str(t) +'\t' + str(Nx)  + '\t' +str(Np) + '\n')


#print np.shape(X),np.shape(z)
#print np.shape(latm.atm_data.a[0]),np.shape(latm.atm_data.a[1])

##########################################################################	
#bomb_flag=0
#random_bomb_place(bomb_num,ep_bomb,Nx,Np,bomb_flag)
#bomb_flag=1
###########################################################################


if (SAVEFIG==1 or PLOTFIG!=0):
	if (PLOT_TYPE==1):
		plot_contour(1,latm.atm_data.a[0],latm.atm_data.a[1],t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Solution T","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)	
		plot_contour(2,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Solution u","Solutuion $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)
		plot_contour(3,latm.atm_data.a[0],latm.atm_data.rqs,t,latm.atm_data.nt,X,P,X,P,[],[],"Solution q","Noise $q_s$",FIG_FILES+'_qs',SAVEFIG,PLOTFIG)
		#plot_contour(3,latm.atm_data.phi,latm.atm_data.phix,t,latm.atm_data.nt,X,P,X[1:Nx+1,1:Np+1],P[1:Nx+1,1:Np+1],colors1,colors2,"Solution $\phi$","Solution $\phi_x$",FIG_FILES+'_phi',SAVEFIG,PLOTFIG)
		#plot_2d(4,latm.atm_data.lda,latm.atm_data.ldax,t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"$\lambda$","$\lambda_x$",FIG_FILES+'_lda',SAVEFIG,PLOTFIG)
		#plot_2d(5,latm.atm_data.ld,CompCond(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"$1/(pA-pB)*d/dx \int_{p_A}^{p_B}u_{tilde} dp$","$1/(pA-pB)*d/dx \int_{p_A}^{p_B}u dp$",FIG_FILES+'_du',SAVEFIG,PLOTFIG)
		#plot_quiver(6,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,"Quiver",FIG_FILES+'_quiver',SAVEFIG,PLOTFIG)
		if (SOURCE==1):
			plot_contour(6,np.abs(latm.atm_data.a[0]-Texac),np.abs(latm.atm_data.a[1]-qexac),t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Error T","Error q",FIG_FILES+'_error',SAVEFIG,PLOTFIG)
			plot_contour(7,Texac,qexac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Exact T","Exact q",FIG_FILES+'_exact',SAVEFIG,PLOTFIG)
			plot_contour(8,abs(latm.atm_data.a[2]-uexac),abs(latm.atm_data.w-wexac),t,latm.atm_data.nt,X,P,X,P,[],[],"Error u","Error $\omega$",FIG_FILES+'_(u,w)_error',SAVEFIG,PLOTFIG)	
			plot_contour(9,uexac,wexac,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Exact u","Exact $\omega$",FIG_FILES+'_(u,w)_exact',SAVEFIG,PLOTFIG)	
	elif (PLOT_TYPE==2):
		plot_contour(1,latm.atm_data.a[0,:,Np/2:Np+2],latm.atm_data.a[1,:,Np/2:Np+2]-q_s(X[:,Np/2:Np+2],P[:,Np/2:Np+2]),t,latm.atm_data.nt,X[:,Np/2:Np+2],z[:,Np/2:Np+2],X[:,Np/2:Np+2],z[:,Np/2:Np+2],colors5,colors2b,"Solution T","Solution q-q_s",FIG_FILES+'_zoom',SAVEFIG,PLOTFIG)	
		plot_contour(21,latm.atm_data.a[0,:,Np/2:Np+2],latm.atm_data.a[1,:,Np/2:Np+2],t,latm.atm_data.nt,X[:,Np/2:Np+2],z[:,Np/2:Np+2],X[:,Np/2:Np+2],z[:,Np/2:Np+2],colors5,colors2,"Solution T","Solution q",FIG_FILES+'_zoom2',SAVEFIG,PLOTFIG)	
		plot_contour(2,latm.atm_data.a[0],latm.atm_data.a[1],t,latm.atm_data.nt,X,z,X,z,colors7,colors2,"Solution T","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)	
		plot_contour(3,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Solutuion $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)
		plot_contour(4,latm.atm_data.a[1]-q_s(X,P),q_s(X,P),t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution diff_q-q_s","Solutuion q_s",FIG_FILES+'_qs',SAVEFIG,PLOTFIG)
		plot_contour_single(5,latm.atm_data.phix,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi'_x$",FIG_FILES+'_phi_x_prime',SAVEFIG,PLOTFIG)
		#plot_2d(4,latm.atm_data.ld,CompCond(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"Before the projection","After the projection",FIG_FILES+'_du',SAVEFIG,PLOTFIG)
		plot_2d(6,latm.atm_data.a[0,:,Np-1],latm.atm_data.a[1,:,Np-1],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant1',SAVEFIG,PLOTFIG)
		plot_2d(6,latm.atm_data.a[0,:,Np-2],latm.atm_data.a[1,:,Np-2],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant2',SAVEFIG,PLOTFIG)
		plot_2d(6,latm.atm_data.a[0,:,Np-3],latm.atm_data.a[1,:,Np-3],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant3',SAVEFIG,PLOTFIG)
		plot_2d(7,latm.atm_data.a[0,:,Np-4],latm.atm_data.a[1,:,Np-4],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant4',SAVEFIG,PLOTFIG)
		plot_2d(7,latm.atm_data.a[0,:,Np-5],latm.atm_data.a[1,:,Np-5],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant5',SAVEFIG,PLOTFIG)
		plot_2d(7,latm.atm_data.a[0,:,Np-6],latm.atm_data.a[1,:,Np-6],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant6',SAVEFIG,PLOTFIG)
		plot_2d(7,latm.atm_data.a[0,:,Np-7],latm.atm_data.a[1,:,Np-7],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant7',SAVEFIG,PLOTFIG)
		plot_2d(7,latm.atm_data.a[0,:,Np-8],latm.atm_data.a[1,:,Np-8],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant8',SAVEFIG,PLOTFIG)
		plot_2d(7,latm.atm_data.a[0,:,Np-9],latm.atm_data.a[1,:,Np-9],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant9',SAVEFIG,PLOTFIG)
		#plot_2d(5,latm.atm_data.a[2,:,Np-10],latm.atm_data.w[:,Np-10],t,latm.atm_data.nt,X[:,0],X[:,0],"u at z="+str(int(z[1,Np-15])),"w at z="+str(int(z[1,Np-10])),FIG_FILES+'_u_at_p_constant1',SAVEFIG,PLOTFIG)
		#plot_2d(5,latm.atm_data.a[2,:,Np-15],latm.atm_data.w[:,Np-15],t,latm.atm_data.nt,X[:,0],X[:,0],"u at z="+str(int(z[1,Np-15])),"w at z="+str(int(z[1,Np-15])),FIG_FILES+'_u_at_p_constant2',SAVEFIG,PLOTFIG)
		#plot_2d(5,latm.atm_data.a[2,:,Np-5],latm.atm_data.w[:,Np-5],t,latm.atm_data.nt,X[:,0],X[:,0],"u at z="+str(int(z[1,Np-15])),"w at z="+str(int(z[1,Np-5])),FIG_FILES+'_u_at_p_constant3',SAVEFIG,PLOTFIG)
		#plot_2d(6,latm.atm_data.a[0,Nx/4,:],latm.atm_data.a[1,Nx/4,:],t,latm.atm_data.nt,P[Nx/4,:],P[Nx/4,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_1',SAVEFIG,PLOTFIG)
		#plot_2d(7,latm.atm_data.a[0,Nx/2-10,:],latm.atm_data.a[1,Nx/2-10,:],t,latm.atm_data.nt,P[Nx/2-10,:],P[Nx/2-10,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_2',SAVEFIG,PLOTFIG)
		#plot_2d(8,latm.atm_data.a[0,Nx/2+10,:],latm.atm_data.a[1,Nx/2+10,:],t,latm.atm_data.nt,P[Nx/2+10,:],P[Nx/2+10,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_4',SAVEFIG,PLOTFIG)
		#plot_2d(9,latm.atm_data.a[0,3*Nx/4,:],latm.atm_data.a[1,3*Nx/4,:],t,latm.atm_data.nt,P[3*Nx/4,:],P[3*Nx/4,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_5',SAVEFIG,PLOTFIG)
		#plot_2d(10,latm.atm_data.a[0,Nx/2,:],latm.atm_data.a[1,Nx/2,:],t,latm.atm_data.nt,P[Nx/2,:],P[Nx/2,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_3',SAVEFIG,PLOTFIG)		
		#plot_contour_single(1,latm.atm_data.pm,t,latm.atm_data.nt,X,z,[],"Pressure",FIG_FILES+'pressure',SAVEFIG,PLOTFIG)
#		plot_curve(6,latm.atm_data.a[0],latm.atm_data.a[1],t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Exact solution q",FIG_FILES+'_T_level',SAVEFIG,PLOTFIG)
#		plot_curve(7,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Exact solution $\omega$",FIG_FILES+'_u_level',SAVEFIG,PLOTFIG)								
#		plot_2d(8,CompCond2(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),CompCond(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"$\int_{p_A}^{p_B}u dp$","$1/(pB-pA) d/dx \int_{p_A}^{p_B}u dp$",FIG_FILES+'_du2',SAVEFIG,PLOTFIG)
#		plot_contour(9,latm.atm_data.a[1],latm.atm_data.a[1]-q_s(X,P),t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution q","q - q_s",FIG_FILES+'_sat',SAVEFIG,PLOTFIG)
#		plot_contour(10,latm.atm_data.a[0] ,300.0 - ((1-P/1000.0))*50.0,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"T","T_bar",FIG_FILES+'_rel0',SAVEFIG,PLOTFIG)				
#		errT=((latm.atm_data.a[0,1:Nx+1,1:Np+1]-(300.0 - (1-P[1:Nx+1,1:Np+1]/1000.0)*50.0))**2)*latm.atm_data.ac
#		errT=np.sqrt(np.sum(errT))/np.sqrt(np.sum(latm.atm_data.a[0,1:Nx+1,1:Np+1]**2.0*latm.atm_data.ac))
#		print errT
		if (SOURCE==1):
			plot_curve(11,latm.atm_data.a[0],np.abs(latm.atm_data.a[0]-Texac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Error T",FIG_FILES+'_Terror',SAVEFIG,PLOTFIG)
			plot_curve(12,latm.atm_data.a[2],np.abs(latm.atm_data.a[2]-uexac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution u","Error u",FIG_FILES+'_Uerror',SAVEFIG,PLOTFIG)
			plot_curve(13,latm.atm_data.w,np.abs(latm.atm_data.w-wexac),t,latm.atm_data.nt,X,z,X,z,[],[],"Solution $\omega$","Error $\omega$",FIG_FILES+'_werror',SAVEFIG,PLOTFIG)	
#			plot_contour(6,np.abs(latm.atm_data.a[0]-Texac),np.abs(latm.atm_data.a[1]-qexac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Solution q",FIG_FILES+'_error',SAVEFIG,PLOTFIG)
#			plot_contour(7,Texac,qexac,t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Solution q",FIG_FILES+'_exact',SAVEFIG,PLOTFIG)
#			plot_contour(8,abs(latm.atm_data.a[2]-uexac),abs(latm.atm_data.w-wexac),t,latm.atm_data.nt,X,P,X,P,[],[],"Error u","Error $\omega$",FIG_FILES+'_(u,w)_error',SAVEFIG,PLOTFIG)	
#			plot_contour(9,uexac,wexac,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Exact u","Exact $\omega$",FIG_FILES+'_(u,w)_exact',SAVEFIG,PLOTFIG)
	elif (PLOT_TYPE==3):
		plot_periodic(1,x0,xf,Nx,Np,latm.atm_data.a[0],latm.atm_data.a[1],t,latm.atm_data.nt,X,z,colors1,colors2,"Solution T","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)
		plot_periodic(2,x0,xf,Nx,Np,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,colors3,colors4,"Solution u","Solution $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)
	elif (PLOT_TYPE==4):
		plot_contour(1,latm.atm_data.a[0]-T_init(X,P),latm.atm_data.a[1],t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T'","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)
		plot_contour(2,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Solution $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)
	elif (PLOT_TYPE==5):
		plot_contour_single(1,latm.atm_data.phi,t,latm.atm_data.nt,X,z,colors1,"Solution $\phi'$",FIG_FILES+'_phi_prime',SAVEFIG,PLOTFIG)
		plot_contour_single(2,latm.atm_data.phix,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi'_x$",FIG_FILES+'_phi_x_prime',SAVEFIG,PLOTFIG)
		plot_contour_single(3,latm.atm_data.phitemp,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi$",FIG_FILES+'_phi',SAVEFIG,PLOTFIG)
		plot_contour_single(4,latm.atm_data.phixtemp,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi_x$",FIG_FILES+'_phi_x',SAVEFIG,PLOTFIG)
		plot_contour_single(7,latm.atm_data.a[2],t,latm.atm_data.nt,X,z,colors1,"Solution $u$",FIG_FILES+'_u',SAVEFIG,PLOTFIG)
		plot_2d_single(5,latm.atm_data.lda,t,latm.atm_data.nt,X[1:Nx+1,0],"bar $\phi$",FIG_FILES+'_bar_phi',SAVEFIG,PLOTFIG)
		plot_2d_single(6,latm.atm_data.ldax,t,latm.atm_data.nt,X[1:Nx+1,0],"bar $\phi_x$",FIG_FILES+'_bar_phi_x',SAVEFIG,PLOTFIG)
		if (SOURCE==1):
			plot_curve(1,latm.atm_data.a[0],Texac,t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Exact solution T",FIG_FILES+'_T',SAVEFIG,PLOTFIG)
			plot_curve(2,latm.atm_data.a[2],uexac,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Exact solution u",FIG_FILES+'_u',SAVEFIG,PLOTFIG)
			plot_curve(3,latm.atm_data.w,wexac,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution $\omega$","Exact solution $\omega$",FIG_FILES+'_w',SAVEFIG,PLOTFIG)
	elif (PLOT_TYPE==6):
		plot_contour(1,latm.atm_data.a[0],Texac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Solution T","Exact solution T",FIG_FILES+'_T',SAVEFIG,PLOTFIG)
		plot_contour(2,latm.atm_data.a[1],qexac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Solution q","Exact solution q",FIG_FILES+'_q',SAVEFIG,PLOTFIG)
		plot_contour(3,latm.atm_data.a[2],uexac,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Solution u","Exact solution u",FIG_FILES+'_u',SAVEFIG,PLOTFIG)
		plot_contour(4,latm.atm_data.w,wexac,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Solution $\omega$","Exact solution $\omega$",FIG_FILES+'_w',SAVEFIG,PLOTFIG)

	
#for i in range(Nx):
	#for j in range(Np):
		#latm.atm_data.a[:,i+1,j+1]=latm.atm_data.a[:,i+1,j+1]+sin(X[i,j]*pi/xf)

################################################################################################
###result 1 Average with count = 5 and Nx,Np = 100
latm.atm_data.a[0,1:Nx+1,:] = (latm.atm_data.a[0,0:Nx,:] + latm.atm_data.a[0,1:Nx+1,:])/2.0
latm.atm_data.a[2,1:Nx+1,:] = (latm.atm_data.a[2,0:Nx,:] + latm.atm_data.a[2,1:Nx+1,:])/2.0

latm.atm_data.w[1:Nx+1,:] = (latm.atm_data.w[0:Nx,:] + latm.atm_data.w[1:Nx+1,:])/2.0
latm.atm_data.phix[1:Nx,:] = (latm.atm_data.phix[0:Nx-1,:] + latm.atm_data.phix[1:Nx,:])/2.0	
latm.atm_data.phix[Nx-1,:] = (latm.atm_data.phix[Nx-1,:] + latm.atm_data.phix[Nx-2,:] + latm.atm_data.phix[Nx-3,:] + latm.atm_data.phix[Nx-4,:])/4.0
latm.atm_data.phix[Nx-2,:] = (latm.atm_data.phix[Nx-5,:] + latm.atm_data.phix[Nx-2,:] + latm.atm_data.phix[Nx-3,:] + latm.atm_data.phix[Nx-4,:])/4.0
################################################################################################

latm.fluxes.borders1()
latm.fluxes.borders3()

print "#### Computation..."

#plt.figure(fig)
plt.clf()
plt.plot(X[1:Nx,0],latm.atm_data.a[1,1:Nx,Np-2],'bs', X[1:Nx,0], temp_qs[1:Nx,Np-2], 'g^')
#plt.ylim((0.0,0.015 ))
plt.legend(['q_initial','q_s'])
plt.savefig('aaa.png')


errT = 0.0
count = 0

count1 = energy_count-1
bomb_count=1
bomb_flag=0
icnt = 0
t_before=time.time()
print tf
#print "X=",X,"dx=",latm.atm_data.dx
#print "P=",P,"dp=",latm.atm_data.dp

if (PLOT_energy==1):
	fbin_T=open('L2_T.dat', 'wb')
	fbin_q=open('L2_q.dat', 'wb')
	fbin_u=open('L2_u.dat', 'wb')
	fbin_w=open('L2_w.dat', 'wb')
	fbin_time=open('L2_time.dat', 'wb')
	fbin_T_profile=open('T_profile.dat', 'wb')
	fbin_q_profile=open('q_profile.dat', 'wb')
	fbin_u_profile=open('u_profile.dat', 'wb')
	fbin_w_profile=open('w_profile.dat', 'wb')
	fbin_T.close()
	fbin_q.close()
	fbin_u.close()
	fbin_w.close()
	fbin_time.close()
	fbin_T_profile.close()
	fbin_q_profile.close()
	fbin_u_profile.close()
	fbin_w_profile.close()
	
	
	
while t<tf-dt/2.0:


	count = count + 1
	count1 = count1 + 1
	
	latm.atm_data.phix[1:Nx,:] = (latm.atm_data.phix[0:Nx-1,:] + latm.atm_data.phix[1:Nx,:])/2.0	
	latm.atm_data.phix[Nx-1,:] = (latm.atm_data.phix[Nx-1,:] + latm.atm_data.phix[Nx-2,:] + latm.atm_data.phix[Nx-3,:] + latm.atm_data.phix[Nx-4,:])/4.0
	latm.atm_data.phix[Nx-2,:] = (latm.atm_data.phix[Nx-5,:] + latm.atm_data.phix[Nx-2,:] + latm.atm_data.phix[Nx-3,:] + latm.atm_data.phix[Nx-4,:])/4.0					
	latm.atm_data.a[2,1:Nx+1,:] = (latm.atm_data.a[2,0:Nx,:] + latm.atm_data.a[2,1:Nx+1,:])/2.0
	latm.atm_data.w[1:Nx+1,:] = (latm.atm_data.w[0:Nx,:] + latm.atm_data.w[1:Nx+1,:])/2.0			
	
	if (PLOT_energy==1 and count1 == energy_count):
		[L2_T,L2_q,L2_u,L2_w] = L2_energy(Nx,Np,latm.atm_data.ac,latm.atm_data.a[0],latm.atm_data.a[1],latm.atm_data.a[2],latm.atm_data.w)
		#latm.atm_data.a[1,1:Nx+1,:] = (latm.atm_data.a[1,0:Nx,:] + latm.atm_data.a[1,1:Nx+1,:])/2.0
		count1 = 0
		bin_T = struct.pack('d', L2_T)
		bin_q = struct.pack('d', L2_q)
		bin_u = struct.pack('d', L2_u)
		bin_w = struct.pack('d', L2_w)
		bin_time = struct.pack('d', t)
		
		fbin_T=open('L2_T.dat', 'ab')
		fbin_q=open('L2_q.dat', 'ab')
		fbin_u=open('L2_u.dat', 'ab')
		fbin_w=open('L2_w.dat', 'ab')
		fbin_time=open('L2_time.dat', 'ab')
		fbin_T_profile=open('T_profile.dat', 'ab')
		fbin_q_profile=open('q_profile.dat', 'ab')
		fbin_u_profile=open('u_profile.dat', 'ab')
		fbin_w_profile=open('w_profile.dat', 'ab')
	
		fbin_T.write(bin_T)
		fbin_q.write(bin_q)
		fbin_u.write(bin_u)
		fbin_w.write(bin_w)
		fbin_time.write(bin_time)
		
		for jk in range(Nx):
			fbin_T_profile.write(latm.atm_data.a[0,jk+1,Np-5])
			fbin_q_profile.write(latm.atm_data.a[1,jk+1,Np-5])
			fbin_u_profile.write(latm.atm_data.a[2,jk+1,Np-5])
			fbin_w_profile.write(latm.atm_data.w[jk+1,Np-5])
		
		fbin_T.close()
		fbin_q.close()
		fbin_u.close()
		fbin_w.close()
		fbin_time.close()
		fbin_T_profile.close()
		fbin_q_profile.close()
		fbin_u_profile.close()
		fbin_w_profile.close()
		
	if (count == AVG1):
		
		#################################################################################################
		# result 1 Average with count = 5 and Nx,Np = 100
		#################################################################################################
		#latm.atm_data.a[:,1:Nx+1,:] = (latm.atm_data.a[:,0:Nx,:] + latm.atm_data.a[:,1:Nx+1,:])/2.0
		
		latm.atm_data.a[0,1:Nx+1,:] = (latm.atm_data.a[0,0:Nx,:] + latm.atm_data.a[0,1:Nx+1,:])/2.0

				
		
		#latm.atm_data.a[0,2:Nx+1,:] = (latm.atm_data.a[0,0:Nx-1,:] + latm.atm_data.a[0,1:Nx,:] + latm.atm_data.a[0,2:Nx+1,:] )/3.0
		#latm.atm_data.a[2,2:Nx+1,:] = (latm.atm_data.a[2,0:Nx-1,:] + latm.atm_data.a[2,1:Nx,:] + latm.atm_data.a[2,2:Nx+1,:] )/3.0

		#latm.atm_data.w[2:Nx+1,:] = (latm.atm_data.w[0:Nx-1,:] + latm.atm_data.w[1:Nx,:] + latm.atm_data.w[2:Nx+1,:] )/3.0
		#latm.atm_data.phix[1:Nx,:] = (latm.atm_data.phix[0:Nx-1,:] + latm.atm_data.phix[1:Nx,:])/2.0	
		#latm.atm_data.phix[Nx-1,:] = (latm.atm_data.phix[Nx-1,:] + latm.atm_data.phix[Nx-2,:] + latm.atm_data.phix[Nx-3,:] + latm.atm_data.phix[Nx-4,:])/4.0
		#latm.atm_data.phix[Nx-2,:] = (latm.atm_data.phix[Nx-5,:] + latm.atm_data.phix[Nx-2,:] + latm.atm_data.phix[Nx-3,:] + latm.atm_data.phix[Nx-4,:])/4.0				
		



		#################################################################################################

		#################################################################################################
		# result 2 Average with count = 8 and Nx,Np = 200
		#################################################################################################
		#latm.atm_data.a[0,1:2,:] = (latm.atm_data.a[0,0:1,:] + latm.atm_data.a[0,1:2,:])/2.0
		#latm.atm_data.a[0,3:Nx+1,:] = (latm.atm_data.a[0,3:Nx+1,:] + latm.atm_data.a[0,2:Nx,:] + latm.atm_data.a[0,1:Nx-1,:] + latm.atm_data.a[0,0:Nx-2,:])/4.0

		#latm.atm_data.a[2,1:2,:] = (latm.atm_data.a[2,0:1,:] + latm.atm_data.a[2,1:2,:])/2.0
		#latm.atm_data.a[2,3:Nx+1,:] = (latm.atm_data.a[2,3:Nx+1,:] + latm.atm_data.a[2,2:Nx,:] + latm.atm_data.a[2,1:Nx-1,:] + latm.atm_data.a[2,0:Nx-2,:])/4.0

		#latm.atm_data.w[1:2,:] = (latm.atm_data.w[0:1,:] + latm.atm_data.w[1:2,:])/2.0
		#latm.atm_data.w[3:Nx+1,:] = (latm.atm_data.w[3:Nx+1,:] + latm.atm_data.w[2:Nx,:] + latm.atm_data.w[1:Nx-1,:] + latm.atm_data.w[0:Nx-2,:])/4.0

		#latm.atm_data.phix[1:2,:] = (latm.atm_data.phix[0:1,:] + latm.atm_data.phix[1:2,:])/2.0
		#latm.atm_data.phix[3:Nx,:] = (latm.atm_data.phix[3:Nx,:] + latm.atm_data.phix[2:Nx-1,:] + latm.atm_data.phix[1:Nx-2,:] + latm.atm_data.phix[0:Nx-3,:])/4.0
		#################################################################################################		
		
		latm.fluxes.borders1()
		latm.fluxes.borders3()
		count = 0
		print count
	
	#for i in range(Nx+2):
	#	for j in range(Np+2):
	#		latm.atm_data.a[2]_2[i,j]=u_ex(X[i,j],P[i,j],t)
	if (SOURCE==1):
		latm.atm_data.s=np.copy(latm.atm_data.sdt)
		latm.atm_data.sdt[0]=ST(X,P,t+dt)
		latm.atm_data.sdt2[0]=ST(X,P,t+dt/2.0)
		latm.atm_data.sdt[1]=Sq(X,P,t+dt)
		latm.atm_data.sdt2[1]=Sq(X,P,t+dt/2.0)
		latm.atm_data.sdt[2]=Su(X,P,t+dt)	
		latm.atm_data.sdt2[2]=Su(X,P,t+dt/2.0)
   
	if (time_method==1):
		print "Euler1..."
		Euler1()
	elif (time_method==2):
		print "Euler2..."
		Euler2()
	elif (time_method==4):
		RK4Source()
	elif (time_method==5):
		#print "Runge Kutta order4, 2..."
		Rk4_2()
	elif (time_method==6):
		#print "Runge Kutta order4, 2..."
		#latm.atm_data.w=w_ex(X,P,t)
		latm.time_method.rk4_3()
		#latm.atm_data.w=w_ex(X,P,t+dt)


	random_qs(ep,Nx,Np) 
	
	
	
	########################################################################################
	
	if (bombing == 1):
		random_bomb_place(bomb_num,ep_bomb,Nx,Np,bomb_flag)
	#bomb_count=bomb_count+1
	
	########################################################################################
	
	
	
	
	t=t+dt
	latm.atm_data.nt=latm.atm_data.nt+1
	print "t=",t
    
	if (SOURCE==1):
		# TEST WITHTOUT U and w
		####################################
		#latm.atm_data.a[2]=u_ex(X2,P2,t)
		#latm.atm_data.w=w_ex(X3,P3,t)
		#latm.atm_data.a[0]=T_ex(X,P,t)
		#latm.atm_data.a[1]=q_ex(X,P,t)
		####################################
		Texac=T_ex(X,P,t)
		qexac=q_ex(X,P,t)
		uexac=u_ex(X,P,t)
		wexac=w_ex(X,P,t)

		#print wexac, latm.atm_data.w,X3, P3,latm.atm_data.ac3
		
		#print latm.atm_data.w
		
		[errT,errq,erru,errw]=L2_norm(Nx,Np,latm.atm_data.ac,latm.atm_data.a[0],latm.atm_data.a[1],latm.atm_data.a[2],latm.atm_data.w,Texac,qexac,uexac,wexac)
		
		print 't=', t, " log(dx)=",np.log10(latm.atm_data.dx)," log(dp)=",np.log10(max(latm.atm_data.dp)), \
		" log(errT)=",np.log10(errT)," log(errq)=",np.log10(errq)," log(erru)=",np.log10(erru)," log(errw)=",np.log10(errw)

#		print 'let do this',latm.atm_data.a[0],latm.atm_data.a[2],latm.atm_data.w
		

#		print " errT=",errT," errq=",errq," erru=",erru," errw=",errw
#	else:
		#print 't=', t
		#print 'max A',  np.max(latm.atm_data.a)
		#print 'max w', np.max(latm.atm_data.w)
		
		#print "T=",latm.atm_data.a[0]
		#print "q=",latm.atm_data.a[1]
		#print "u=",latm.atm_data.a[2]
		#print "w=",latm.atm_data.w
		#print "Look at me Dude!!",latm.atm_data.s
	if (DEBUG ==1) :
		plt.figure(1)
		plt.clf()
		plt.plot(X[1:Nx+1,1],latm.atm_data.lda)
		plt.title('lmbda at t='+str(t))
		plt.savefig('lda_'+str(t)+'.png')
		plt.figure(1)
		plt.clf()
		plt.plot(X[1:Nx+1,1],latm.atm_data.ldax)
		plt.title('lmbdax at t='+str(t))
		plt.savefig('ldax_'+str(t)+'.png')
		plt.figure(1)
		plt.clf()
		plt.quiver(X,P,latm.atm_data.a[2],latm.atm_data.w)
		plt.title('quiver at t='+str(t))
		plt.savefig('quiver_'+str(t)+'.png')
		plt.figure(1)
		plt.clf()
		div=np.zeros((Nx,Np),'d',order='f')
		for i in range(Nx):
			for j in range(Np):
				div[i,j]=(latm.atm_data.a[2,i+2,j+1]-latm.atm_data.a[2,i,j+1])/(2*latm.atm_data.dx)+(latm.atm_data.w[i+1,j+2]-latm.atm_data.w[i+1,j])/(latm.atm_data.dp[i]+latm.atm_data.dp[i]+1)
		plt.contourf(X[1:Nx+1,1:Np+1],P[1:Nx+1,1:Np+1],div)
		plt.colorbar(extend="both",format="%.2e")
		plt.title('divergence at t='+str(t))
		plt.savefig('divergence_'+str(t)+'.png')

	if (PLOT_ERROR==1 and latm.atm_data.nt%PERIODE_ERROR==0):
		errT_t[ind_t]=errT/np.sqrt(latm.atm_data.dp*np.inner(Texac,Texac))
		errq_t[ind_t]=errq/np.sqrt(latm.atm_data.dp*np.inner(qexac,qexac))
		erru_t[ind_t]=erru/np.sqrt(latm.atm_data.dp*np.inner(Texac,Texac))
		errw1_t[ind_t]=errq/np.sqrt(latm.atm_data.dp*np.inner(qexac,qexac))
		errw2_t[ind_t]=errq/np.sqrt(latm.atm_data.dp*np.inner(qexac,qexac))
		t_t[ind_t]=t
		ind_t=ind_t+1


	
	if (SAVEFIG==1 or PLOTFIG!=0):
		if (latm.atm_data.nt%PERIODE_FIG==0):
			#FILE_ID.write(str(errT) + '\n')
			if (PLOT_TYPE==1):
				plot_contour(1,latm.atm_data.a[0],latm.atm_data.a[1],t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Solution T","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)	
				plot_contour(2,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Solution u","Solutuion $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)
				plot_contour(3,latm.atm_data.a[0],latm.atm_data.rqs,t,latm.atm_data.nt,X,P,X,P,[],[],"Solution q","Noise $q_s$",FIG_FILES+'_qs',SAVEFIG,PLOTFIG)
				#plot_contour(3,latm.atm_data.phi,latm.atm_data.phix,t,latm.atm_data.nt,X,P,X[1:Nx+1,1:Np+1],P[1:Nx+1,1:Np+1],colors1,colors2,"Solution $\phi$","Solution $\phi_x$",FIG_FILES+'_phi',SAVEFIG,PLOTFIG)
				#plot_2d(4,latm.atm_data.lda,latm.atm_data.ldax,t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"$\lambda$","$\lambda_x$",FIG_FILES+'_lda',SAVEFIG,PLOTFIG)
				#plot_2d(5,latm.atm_data.ld,CompCond(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"$d/dx \int_{p_A}^{p_B}u_{tilde} dp$","$d/dx \int_{p_A}^{p_B}u dp$",FIG_FILES+'_du',SAVEFIG,PLOTFIG)
				#plot_quiver(6,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,"Quiver",FIG_FILES+'_quiver',SAVEFIG,PLOTFIG)
				if (SOURCE==1):
					plot_contour(6,np.abs(latm.atm_data.a[0]-Texac),np.abs(latm.atm_data.a[1]-qexac),t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Error T","Error q",FIG_FILES+'_error',SAVEFIG,PLOTFIG)
					plot_contour(7,Texac,qexac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Exact T","Exact q",FIG_FILES+'_exact',SAVEFIG,PLOTFIG)
					plot_contour(8,abs(latm.atm_data.a[2]-uexac),abs(latm.atm_data.w-wexac),t,latm.atm_data.nt,X,P,X,P,[],[],"Error u","Error $\omega$",FIG_FILES+'_(u,w)_error',SAVEFIG,PLOTFIG)	
					plot_contour(9,uexac,wexac,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Exact u","Exact $\omega$",FIG_FILES+'_(u,w)_exact',SAVEFIG,PLOTFIG)		
			elif (PLOT_TYPE==2):
				plot_contour(1,latm.atm_data.a[0,:,Np/2:Np+2],latm.atm_data.a[1,:,Np/2:Np+2]-q_s(X[:,Np/2:Np+2],P[:,Np/2:Np+2]),t,latm.atm_data.nt,X[:,Np/2:Np+2],z[:,Np/2:Np+2],X[:,Np/2:Np+2],z[:,Np/2:Np+2],colors5,colors2b,"Solution T","Solution q-q_s",FIG_FILES+'_zoom',SAVEFIG,PLOTFIG)	
				plot_contour(1,latm.atm_data.a[0,:,Np/2:Np+2],latm.atm_data.a[1,:,Np/2:Np+2],t,latm.atm_data.nt,X[:,Np/2:Np+2],z[:,Np/2:Np+2],X[:,Np/2:Np+2],z[:,Np/2:Np+2],colors5,colors2,"Solution T","Solution q",FIG_FILES+'_zoom2',SAVEFIG,PLOTFIG)	
				#plot_contour(1,latm.atm_data.a[0,:,3*Np/4:Np+2],latm.atm_data.a[1,:,3*Np/4:Np+2],t,latm.atm_data.nt,X[:,3*Np/4:Np+2],z[:,3*Np/4:Np+2],X[:,3*Np/4:Np+2],z[:,3*Np/4:Np+2],colors5,colors2,"Solution T","Solution q",FIG_FILES+'_zoom',SAVEFIG,PLOTFIG)	
				plot_contour(2,latm.atm_data.a[0],latm.atm_data.a[1],t,latm.atm_data.nt,X,z,X,z,colors7,colors2,"Solution T","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)	
				plot_contour(3,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Solutuion $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)
				plot_contour(4,latm.atm_data.a[1]-q_s(X,P),q_s(X,P),t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution diff_q-q_s","Solutuion q_s",FIG_FILES+'_qs',SAVEFIG,PLOTFIG)
				plot_contour_single(5,latm.atm_data.phix,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi'_x$",FIG_FILES+'_phi_x_prime',SAVEFIG,PLOTFIG)
				plot_2d(6,latm.atm_data.a[0,:,Np-1],latm.atm_data.a[1,:,Np-1],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant1',SAVEFIG,PLOTFIG)
				plot_2d(6,latm.atm_data.a[0,:,Np-2],latm.atm_data.a[1,:,Np-2],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant2',SAVEFIG,PLOTFIG)
				plot_2d(6,latm.atm_data.a[0,:,Np-3],latm.atm_data.a[1,:,Np-3],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant3',SAVEFIG,PLOTFIG)
				plot_2d(7,latm.atm_data.a[0,:,Np-4],latm.atm_data.a[1,:,Np-4],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant4',SAVEFIG,PLOTFIG)
				plot_2d(7,latm.atm_data.a[0,:,Np-5],latm.atm_data.a[1,:,Np-5],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant5',SAVEFIG,PLOTFIG)
				plot_2d(7,latm.atm_data.a[0,:,Np-6],latm.atm_data.a[1,:,Np-6],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant6',SAVEFIG,PLOTFIG)
				plot_2d(7,latm.atm_data.a[0,:,Np-7],latm.atm_data.a[1,:,Np-7],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant7',SAVEFIG,PLOTFIG)
				plot_2d(7,latm.atm_data.a[0,:,Np-8],latm.atm_data.a[1,:,Np-8],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant8',SAVEFIG,PLOTFIG)
				plot_2d(7,latm.atm_data.a[0,:,Np-9],latm.atm_data.a[1,:,Np-9],t,latm.atm_data.nt,X[:,0],X[:,0],"T near the mounatain","q near the mounatain",FIG_FILES+'_T_at_p_constant9',SAVEFIG,PLOTFIG)
				#plot_2d(5,latm.atm_data.a[2,:,Np-10],latm.atm_data.w[:,Np-10],t,latm.atm_data.nt,X[:,0],X[:,0],"u at z="+str(int(z[1,Np-15])),"w at z="+str(int(z[1,Np-10])),FIG_FILES+'_u_at_p_constant1',SAVEFIG,PLOTFIG)
				#plot_2d(5,latm.atm_data.a[2,:,Np-15],latm.atm_data.w[:,Np-15],t,latm.atm_data.nt,X[:,0],X[:,0],"u at z="+str(int(z[1,Np-15])),"w at z="+str(int(z[1,Np-15])),FIG_FILES+'_u_at_p_constant2',SAVEFIG,PLOTFIG)
				#plot_2d(5,latm.atm_data.a[2,:,Np-5],latm.atm_data.w[:,Np-5],t,latm.atm_data.nt,X[:,0],X[:,0],"u at z="+str(int(z[1,Np-15])),"w at z="+str(int(z[1,Np-5])),FIG_FILES+'_u_at_p_constant3',SAVEFIG,PLOTFIG)
				########################
				#plot_2d(4,latm.atm_data.ld,CompCond(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"Before the projection","After the projection",FIG_FILES+'_du',SAVEFIG,PLOTFIG)
				#plot_2d(5,latm.atm_data.a[0,:,Np],latm.atm_data.a[1,:,Np-2],t,latm.atm_data.nt,X[:,0],X[:,0],"T_near_mounatain","q_near_mounatain",FIG_FILES+'_T_at_p_constant',SAVEFIG,PLOTFIG)
				#plot_2d(6,latm.atm_data.a[0,Nx/4,:],latm.atm_data.a[1,Nx/4,:],t,latm.atm_data.nt,P[Nx/4,:],P[Nx/4,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_1',SAVEFIG,PLOTFIG)
				#plot_2d(7,latm.atm_data.a[0,Nx/2-10,:],latm.atm_data.a[1,Nx/2-10,:],t,latm.atm_data.nt,P[Nx/2-10,:],P[Nx/2-10,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_2',SAVEFIG,PLOTFIG)
				#plot_2d(8,latm.atm_data.a[0,Nx/2+10,:],latm.atm_data.a[1,Nx/2+10,:],t,latm.atm_data.nt,P[Nx/2+10,:],P[Nx/2+10,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_4',SAVEFIG,PLOTFIG)
				#plot_2d(9,latm.atm_data.a[0,3*Nx/4,:],latm.atm_data.a[1,3*Nx/4,:],t,latm.atm_data.nt,P[3*Nx/4,:],P[3*Nx/4,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_5',SAVEFIG,PLOTFIG)
				#plot_2d(10,latm.atm_data.a[0,Nx/2,:],latm.atm_data.a[1,Nx/2,:],t,latm.atm_data.nt,P[Nx/2,:],P[Nx/2,:],"T_vertical","q_vertical",FIG_FILES+'_T_q_at_x_constant_3',SAVEFIG,PLOTFIG)		
				########################
				plot_curve(8,latm.atm_data.a[0,:,50:],latm.atm_data.a[1,:,50:],t,latm.atm_data.nt,X[:,50:],z[:,50:],X[:,50:],z[:,50:],[],[],"Solution T","Exact solution q",FIG_FILES+'_T_level',SAVEFIG,PLOTFIG)
				plot_curve(9,latm.atm_data.a[2,:,50:],latm.atm_data.w[:,50:],t,latm.atm_data.nt,X[:,50:],z[:,50:],X[:,50:],z[:,50:],np.linspace(5.0,30.0,30),np.linspace(-10.0,10.0,30),"Solution u","Exact solution $\omega$",FIG_FILES+'_u_level',SAVEFIG,PLOTFIG)								
#				plot_2d(8,CompCond2(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),CompCond(Nx,Np,latm.atm_data.dp,latm.atm_data.dx,P,latm.atm_data.a[2]),t,latm.atm_data.nt,X[1:Nx+1,0],X[1:Nx+1,0],"$\int_{p_A}^{p_B}u dp$","$1/(pB-pA) d/dx \int_{p_A}^{p_B}u dp$",FIG_FILES+'_du2',SAVEFIG,PLOTFIG)
#				plot_contour(9,latm.atm_data.a[1],latm.atm_data.a[1]-q_s(X,P),t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution q","q - q_s",FIG_FILES+'_sat',SAVEFIG,PLOTFIG)
#				plot_contour(10,latm.atm_data.a[0] ,300.0 - ((1-P/1000.0))*50.0,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"T","T_bar",FIG_FILES+'_rel0',SAVEFIG,PLOTFIG)				
#				errT=((latm.atm_data.a[0,1:Nx+1,1:Np+1]-(300.0 - (1-P[1:Nx+1,1:Np+1]/1000.0)*50.0))**2)*latm.atm_data.ac
#				errT=np.sqrt(np.sum(errT))/np.sqrt(np.sum(latm.atm_data.a[0,1:Nx+1,1:Np+1]**2.0*latm.atm_data.ac))
#				print errT
				if (SOURCE==1):
					plot_curve(11,latm.atm_data.a[0],np.abs(latm.atm_data.a[0]-Texac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Error T",FIG_FILES+'_Terror',SAVEFIG,PLOTFIG)
					plot_curve(12,latm.atm_data.a[2],np.abs(latm.atm_data.a[2]-uexac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution u","Error u",FIG_FILES+'_Uerror',SAVEFIG,PLOTFIG)
					plot_curve(13,latm.atm_data.w,np.abs(latm.atm_data.w-wexac),t,latm.atm_data.nt,X,z,X,z,[],[],"Solution $\omega$","Error $\omega$",FIG_FILES+'_werror',SAVEFIG,PLOTFIG)	
#					plot_contour(6,np.abs(latm.atm_data.a[0]-Texac),np.abs(latm.atm_data.a[1]-qexac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Solution q",FIG_FILES+'_error',SAVEFIG,PLOTFIG)
#					plot_contour(7,Texac,qexac,t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Solution q",FIG_FILES+'_exact',SAVEFIG,PLOTFIG)
#					plot_contour(8,abs(latm.atm_data.a[2]-uexac),abs(latm.atm_data.w-wexac),t,latm.atm_data.nt,X,P,X,P,[],[],"Error u","Error $\omega$",FIG_FILES+'_(u,w)_error',SAVEFIG,PLOTFIG)	
#					plot_contour(9,uexac,wexac,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Exact u","Exact $\omega$",FIG_FILES+'_(u,w)_exact',SAVEFIG,PLOTFIG)
			elif (PLOT_TYPE==3):
				plot_periodic(1,x0,xf,Nx,Np,latm.atm_data.a[0],latm.atm_data.a[1],t,latm.atm_data.nt,X,z,colors1,colors2,"Solution T","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)	
				plot_periodic(2,x0,xf,Nx,Np,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,colors3,colors4,"Solution u","Solutuion $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)				
			elif (PLOT_TYPE==4):
				plot_contour(1,latm.atm_data.a[0]-T_bar(X,P),latm.atm_data.a[1],t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T'","Solution q",FIG_FILES,SAVEFIG,PLOTFIG)
				plot_contour(2,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Solution $\omega$",FIG_FILES+'_(u,w)',SAVEFIG,PLOTFIG)	
				plot_curve(3,latm.atm_data.a[0]-T_bar(X,P),latm.atm_data.a[1],t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T'","Exact solution q",FIG_FILES+'_T_level',SAVEFIG,PLOTFIG)
				plot_curve(4,latm.atm_data.a[2],latm.atm_data.w,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Exact solution $omega$",FIG_FILES+'_u_level',SAVEFIG,PLOTFIG)								
			elif (PLOT_TYPE==5):
				plot_contour_single(1,latm.atm_data.phi,t,latm.atm_data.nt,X,z,colors1,"Solution $\phi'$",FIG_FILES+'phi_prime',SAVEFIG,PLOTFIG)
				plot_contour_single(2,latm.atm_data.phix,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi'_x$",FIG_FILES+'_phi_x_prime',SAVEFIG,PLOTFIG)
				plot_contour_single(3,latm.atm_data.phitemp,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi$",FIG_FILES+'_phi',SAVEFIG,PLOTFIG)
				plot_contour_single(4,latm.atm_data.phixtemp,t,latm.atm_data.nt,X[1:Nx+1,1:Np+1],z[1:Nx+1,1:Np+1],colors1,"Solution $\phi_x$",FIG_FILES+'_phi_x',SAVEFIG,PLOTFIG)
				plot_contour_single(7,latm.atm_data.a[2],t,latm.atm_data.nt,X,z,colors1,"Solution $u$",FIG_FILES+'_u',SAVEFIG,PLOTFIG)
				plot_2d_single(5,latm.atm_data.lda,t,latm.atm_data.nt,X[1:Nx+1,0],"bar $\phi$",FIG_FILES+'_bar_phi',SAVEFIG,PLOTFIG)
				plot_2d_single(6,latm.atm_data.ldax,t,latm.atm_data.nt,X[1:Nx+1,0],"bar $\phi_x$",FIG_FILES+'_bar_phi_x',SAVEFIG,PLOTFIG)
				if (SOURCE==1):
					plot_curve(1,latm.atm_data.a[0],Texac,t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"Solution T","Exact solution T",FIG_FILES+'_T',SAVEFIG,PLOTFIG)
					plot_curve(2,latm.atm_data.a[2],uexac,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution u","Exact solution u",FIG_FILES+'_u',SAVEFIG,PLOTFIG)
					plot_curve(3,latm.atm_data.w,wexac,t,latm.atm_data.nt,X,z,X,z,colors3,colors4,"Solution $\omega$","Exact solution $\omega$",FIG_FILES+'_w',SAVEFIG,PLOTFIG)
					plot_curve(4,np.abs(latm.atm_data.a[0]-Texac),np.abs(Texac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"error for T","Exact solution T",FIG_FILES+"error_",SAVEFIG,PLOTFIG)
					#plot_curve(5,np.abs(latm.atm_data.a[2]-uexac),wexac,t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"error for u","$\omega$",FIG_FILES+"error_(u,w)",SAVEFIG,PLOTFIG)
					plot_curve(5,np.abs(latm.atm_data.a[2]-uexac),np.abs(latm.atm_data.w-wexac),t,latm.atm_data.nt,X,z,X,z,colors1,colors2,"error for u","error for $\omega$",FIG_FILES+"error_(u,w)",SAVEFIG,PLOTFIG)
#					plot_curve(6,latm.atm_data.phi,latm.atm_data.phix,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Solution of $\phi$","Solution of $\phi_x$",FIG_FILES+"phi",SAVEFIG,PLOTFIG)
			elif (PLOT_TYPE==6):
				plot_contour(1,latm.atm_data.a[0],Texac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Solution T","Exact solution T",FIG_FILES+'_T',SAVEFIG,PLOTFIG)
				plot_contour(2,latm.atm_data.a[1],qexac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"Solution q","Exact solution q",FIG_FILES+'_q',SAVEFIG,PLOTFIG)
				plot_contour(3,latm.atm_data.a[2],uexac,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Solution u","Exact solution u",FIG_FILES+'_u',SAVEFIG,PLOTFIG)
				plot_contour(4,latm.atm_data.w,wexac,t,latm.atm_data.nt,X,P,X,P,colors3,colors4,"Solution $\omega$","Exact solution $\omega$",FIG_FILES+'_w',SAVEFIG,PLOTFIG)
				plot_contour(5,np.abs(latm.atm_data.a[0]-Texac),Texac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"error for T","Exact solution T",FIG_FILES+"error_",SAVEFIG,PLOTFIG)
				plot_contour(6,np.abs(latm.atm_data.a[1]-qexac),qexac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"error for q","Exact solution q",FIG_FILES+"error_",SAVEFIG,PLOTFIG)
				plot_contour(7,np.abs(latm.atm_data.a[2]-uexac),wexac,t,latm.atm_data.nt,X,P,X,P,colors1,colors2,"error for u","$\omega$",FIG_FILES+"error_(u,w)",SAVEFIG,PLOTFIG)
					
	if (SAVE2TXT_DATA==1):
		if (latm.atm_data.nt%PERIODE_DATA==0):
			save_txt(DATA_DIR,t,Nx,Np,X,P,latm.atm_data.a[0],latm.atm_data.a[1],latm.atm_data.a[2],latm.atm_data.w)


	
t_after=time.time()

t_elapsed=t_after-t_before

print "#### Elapsed time = ",t_elapsed


#if (PLOT_energy==1):
	#fbin_T.close()
	#fbin_q.close()
	#fbin_u.close()
	#fbin_w.close()
	#fbin_time.close()

if (PLOT_ERROR==1):
	plt.figure(1)	
	plt.clf()

	plt.subplot(411)
	plt.plot(t_t,errT_t)
	#plt.xlabel('t')
	plt.title('Error L2 of T over time',fontsize=12)

	plt.subplot(412)
	plt.plot(t_t,errq_t)
	plt.xlabel('t')
	plt.title('Error L2 of q over time',fontsize=12)
	
	plt.subplot(421)
	plt.plot(t_t,errq_t)
	plt.xlabel('t')
	plt.title('Error L2 of u over time',fontsize=12)
	
	plt.subplot(422)
	plt.plot(t_t,errq_t)
	plt.xlabel('t')
	plt.title('Error L2 of w over time',fontsize=12)

	plt.draw()
	plt.savefig('L2_energy.png')
		
if (DEBUG==1):
	print 'You can see the computations in DEBUG/ ...'
if (SAVEFIG==1):
	print 'You can see the graphism in IMG/...'
if (SAVE2TXT_DATA == 1):
	print 'You can see the results in DATA/ ...'

print '######################################'
print '#### End Resolution               ####'
print '######################################'
