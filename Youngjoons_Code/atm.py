# import library
import numpy as np
import libatm as latm
from numpy import arange
#from mlplot import *
import os




# Relative L^2 error T,q,u and \omega
def L2_norm(Nx,Np,ac,T,q,u,w,Texac,qexac,uexac,wexac):
    
    errT=((T[1:Nx+1,1:Np+1]-Texac[1:Nx+1,1:Np+1])**2)*ac
    errq=((q[1:Nx+1,1:Np+1]-qexac[1:Nx+1,1:Np+1])**2)*ac
    erru=((u[1:Nx+1,1:Np+1]-uexac[1:Nx+1,1:Np+1])**2)*ac
    errw=((w[1:Nx+1,1:Np+1]-wexac[1:Nx+1,1:Np+1])**2)*ac
    errT=np.sqrt(np.sum(errT))/np.sqrt(np.sum(Texac[1:Nx+1,1:Np+1]**2.0*ac))
    errq=np.sqrt(np.sum(errq))/np.sqrt(np.sum(qexac[1:Nx+1,1:Np+1]**2.0*ac))
    erru=np.sqrt(np.sum(erru))/np.sqrt(np.sum(uexac[1:Nx+1,1:Np+1]**2.0*ac))
    errw=np.sqrt(np.sum(errw))/np.sqrt(np.sum(wexac[1:Nx+1,1:Np+1]**2.0*ac))
    
    return [errT,errq,erru,errw]

# L^2 norm of T,q,u and \omega
def L2_energy(Nx,Np,ac,T,q,u,w):
    
    L2_T=np.sqrt(np.sum(((T[1:Nx+1,1:Np+1])**2)*ac))
    L2_q=np.sqrt(np.sum(((q[1:Nx+1,1:Np+1])**2)*ac))
    L2_u=np.sqrt(np.sum(((u[1:Nx+1,1:Np+1])**2)*ac))
    L2_w=np.sqrt(np.sum(((w[1:Nx+1,1:Np+1])**2)*ac))
    
    
    return [L2_T,L2_q,L2_u,L2_w]
    
# H^1 norm of T,q,u and \omega
def H1_energy(Nx,Np,ac,T,q,u,w):
    
    gradT=np.zeros((Nx,Np),'d')
    
    
    for i in range(Nx):
		for j in range(Np):
			gradT=(latm.atm_data.m2[i,j,0]*(T[i+2,j+1]-T[i,j+1])+latm.atm_data.m2[i,j,1]*(T[i+1,j+2]-T[i+1,j]))**2 \
						+ ((T[i+1,j+2]-T[i+1,j])/(latm.atm_data.dp[i]+latm.atm_data.dp[i+1]))**2
			gradq=(latm.atm_data.m2[i,j,0]*(q[i+2,j+1]-q[i,j+1])+latm.atm_data.m2[i,j,1]*(q[i+1,j+2]-q[i+1,j]))**2 \
						+ ((q[i+1,j+2]-q[i+1,j])/(latm.atm_data.dp[i]+latm.atm_data.dp[i+1]))**2
			gradu=(latm.atm_data.m2[i,j,0]*(u[i+2,j+1]-u[i,j+1])+latm.atm_data.m2[i,j,1]*(u[i+1,j+2]-u[i+1,j]))**2 \
						+ ((u[i+1,j+2]-u[i+1,j])/(latm.atm_data.dp[i]+latm.atm_data.dp[i+1]))**2
			gradw=(latm.atm_data.m2[i,j,0]*(w[i+2,j+1]-w[i,j+1])+latm.atm_data.m2[i,j,1]*(w[i+1,j+2]-w[i+1,j]))**2 \
						+ ((w[i+1,j+2]-w[i+1,j])/(latm.atm_data.dp[i]+latm.atm_data.dp[i+1]))**2
    H1_T=np.sqrt(np.sum(((T[1:Nx+1,1:Np+1])**2+gradT)*ac))
    H1_q=np.sqrt(np.sum(((q[1:Nx+1,1:Np+1])**2+gradq)*ac))
    H1_u=np.sqrt(np.sum(((u[1:Nx+1,1:Np+1])**2+gradu)*ac))
    H1_w=np.sqrt(np.sum(((w[1:Nx+1,1:Np+1])**2+gradw)*ac))
    
    
    return [H1_T,H1_q,H1_u,H1_w]
"""***************************************************
This is function for meshing of x coordinates
This function send information to libatm, i.e. atm_data.f90 in lflux folder

 Convert the data in python to fortran in the sw_data.f90
 Converted / created data in fortran are :
    sw_data.n
    sw_data.p
    sw_data.pm
    sw_data.dt
    sw_data.dx
    sw_data.theta
    sw_data.debug

Input parameter
    dt : delta t, time step, double
    g :  9.8
    DEBUG = 0 debug mod desactivate , 1 debug mod activate
    x0 : starting point of x
    xf : end point of x
    Nx : the number of disretization for x axis 
    p0 : p_0 in the paper  , 1000
    Np : the number of discretiztion for p axis (similar to y-axis)
    PHYSIC : 
Output parameter
    latm.atm_data.x            - x coordinates for standard mesh, (Nx+1,Np+1) matrix
    latm.atm_data.xm        - x coordinates for T and q, (Nx+2,Np+2) matrix
    latm.atm_data.xm2        - x coordinates for u, (Nx+1,Np+2) matrix
    latm.atm_data.xm3        - x coordinates for w, (Nx+2,Np+1) matrix
    latm.atm_data.p0        - Set input parameter p0
    latm.atm_data.dx        - find delta x using x0 and xf 
    latm.atm_data.dt        - Set input parameter dt
    latm.atm_data.g            - Set input parameter g 
    latm.atm_data.debug        - Set input parameter debug
    latm.atm_data.physics    - Set input parameter physics
    
    If we ratate our matrix 90 degree count-clockwise, then we can match
    our mesh.
***************************************************"""


def data_f90(dt,g,DEBUG,x0,xf,Nx,p0,Np,PHYSIC):

    #    Initialize variables. Then number of discretizations
    latm.atm_data.nx=Nx
    latm.atm_data.np=Np    
    
    #     fucntion zeros is defined in numpy (np)
    #    2 x Nx x Np matrix wih data type 'd' = double precision floating point number
    #    and store option is fortran.
    
    latm.atm_data.fluxhori=np.zeros((3,Nx,Np),'d',order='f')
    latm.atm_data.fluxvert=np.zeros((3,Nx,Np),'d',order='f')
    latm.atm_data.da=np.zeros((3,Nx,Np),'d',order='f')
    latm.atm_data.phi=np.zeros((Nx+2,Np+2),'d',order='f')
    latm.atm_data.qs_temp=np.zeros((Nx+2,Np+2),'d',order='f')
    latm.atm_data.phix=np.zeros((Nx,Np),'d',order='f')
    latm.atm_data.phitemp=np.zeros((Nx,Np),'d',order='f')
    latm.atm_data.phixtemp=np.zeros((Nx,Np),'d',order='f')
    latm.atm_data.x=np.zeros((Nx+1,Np+1),'d',order='f')
	
    latm.atm_data.xm2=np.zeros((Nx+1,Np+2),'d',order='f')
    latm.atm_data.xm3=np.zeros((Nx+2,Np+1),'d',order='f')
    latm.atm_data.vcheck=np.zeros(Nx*(Np+1),'d',order='f')
    latm.atm_data.lda=np.zeros((Nx),'d',order='f')
    latm.atm_data.ldtemp=np.zeros((Nx),'d',order='f')
    latm.atm_data.ldax=np.zeros((Nx),'d',order='f')
    latm.atm_data.ld=np.zeros((Nx),'d',order='f')
    #    get the value p0 and delta x using input parameters
    latm.atm_data.p0=p0
    latm.atm_data.dx = (xf-x0)/float(Nx)

    # for each j, we can give same value of x, xm. That is, x coordinates is same whatever p is.
    for j in range(Np):
        latm.atm_data.x[:,j]=x0+np.arange(Nx+1)*latm.atm_data.dx
        
    latm.atm_data.x[:,Np]=x0+np.arange(Nx+1)*latm.atm_data.dx
    
    # data inside of meshes
    for j in range(Np):
        latm.atm_data.xm2[:,j+1]=(latm.atm_data.x[:,j]+latm.atm_data.x[:,j+1])/2.0
    
    # data on the boundary
    latm.atm_data.xm2[:,0]=latm.atm_data.x[:,0]
    latm.atm_data.xm2[:,Np+1]=latm.atm_data.x[:,Np]
    
    # data inside of meshes
    for i in range(Nx):
        latm.atm_data.xm3[i+1,:]=(latm.atm_data.x[i+1,:]+latm.atm_data.x[i,:])/2.0
    
    # data on the boundary
    latm.atm_data.xm3[0,:]=latm.atm_data.x[0,:]- latm.atm_data.dx/2.0
    latm.atm_data.xm3[Nx+1,:]=latm.atm_data.x[Nx,:] + latm.atm_data.dx/2.0
    
    
    
    latm.atm_data.dt = dt
    latm.atm_data.g = g
    latm.atm_data.debug = DEBUG
    

    latm.atm_data.physic=PHYSIC

"""***************************************************
This is function for meshing of x coordinates
This function send information to libatm, i.e. atm_data.f90 in lflux folder

Input parameter
    pB     -    vector size (Nx+1) containing the value of pB(x) for x.
    pA     -     real number ~ 10, at the top of the atmosphere.
Output parameter
    latm.atm_data.p            - p coordinates for standard mesh, (Nx+1,Np+1) matrix
    latm.atm_data.pm        - p coordinates for T and q, (Nx+2,Np+2) matrix
    latm.atm_data.pm2        - p coordinates for u, (Nx+1,Np+2) matrix
    latm.atm_data.pm3        - p coordinates for w, (Nx+2,Np+1) matrix
    latm.atm_data.dp        - length of each delta p for each corresponding x
                              dim(dp) = Nx+1.
    latm.atm_data.normal    - normal Matrix. Mat = (2,Nx*(Np+1)).
                              1st comp. -> x comp. of normal vercot
                              2st comp. -> y comp. of normal vercot
                              The order is left to right and then 
                              bottom to up. (in the sense of Meshes)
    latm.atm_data.hhori        - length of each dp. It's a vertor with length (Np+1)*Nx
                              This comp. agrees with the normal Matrix 2nd component.
    latm.atm.data.ac        - area for T,q. Mat = (Nx,Np)
    latm.atm_data.ac2        - area for u. Mat = (Nx,Np)
    latm.atm_data.ac3        - area for w. Mat = (Nx,Np+1)
    
    If we ratate our matrix 90 degree count-clockwise, then we can match
    our mesh.
***************************************************"""

def mesh(x0,xf,pB,pA,Nx,Np):
    #initialized
    latm.atm_data.p=np.zeros((Nx+1,Np+1),'d',order='f')
    latm.atm_data.xm=np.zeros((Nx+2,Np+2),'d',order='f')
    latm.atm_data.pm=np.zeros((Nx+2,Np+2),'d',order='f')

    latm.atm_data.pm2=np.zeros((Nx+1,Np+2),'d',order='f')
    latm.atm_data.pm3=np.zeros((Nx+2,Np+1),'d',order='f')
    latm.atm_data.dp=np.zeros((Nx+1),'d',order='f')
    latm.atm_data.normal=np.zeros((2,Nx*(Np+1)),'d',order='f')
    latm.atm_data.hhori=np.zeros(((Np+1)*Nx),'d',order='f')
    latm.atm_data.ac=np.zeros((Nx,Np),'d',order='f')
    latm.atm_data.rqs=np.zeros((Nx+2,Np+2),'d',order='f')

    # center of the mesh
    xc_tri1 = np.zeros((Nx,Np),'d',order='f')
    xc_tri2 = np.zeros((Nx,Np),'d',order='f')
    xc_tri3 = np.zeros((Nx,Np),'d',order='f')
    xc_tri4 = np.zeros((Nx,Np),'d',order='f')
    pc_tri1 = np.zeros((Nx,Np),'d',order='f')
    pc_tri2 = np.zeros((Nx,Np),'d',order='f')
    pc_tri3 = np.zeros((Nx,Np),'d',order='f')
    pc_tri4 = np.zeros((Nx,Np),'d',order='f')
    A_1 = np.zeros((Nx,Np),'d',order='f')
    A_2 = np.zeros((Nx,Np),'d',order='f')
    L = xf-x0
    # standard mesh
    for i in range(Nx+1):
        latm.atm_data.dp[i]=(pB[i]-pA)/Np
        latm.atm_data.p[i,:]=pA+np.arange(Np+1)*latm.atm_data.dp[i]
    
    
    """
    To get the center of mesh
    """
    
    xc_tri1[0:Nx,0:Np] = (latm.atm_data.x[0:Nx  ,0:Np] + latm.atm_data.x[1:Nx+1,0:Np  ] + latm.atm_data.x[0:Nx  ,1:Np+1])/3.0 
    xc_tri2[0:Nx,0:Np] = (latm.atm_data.x[0:Nx  ,0:Np] + latm.atm_data.x[1:Nx+1,0:Np  ] + latm.atm_data.x[1:Nx+1,1:Np+1])/3.0 
    xc_tri3[0:Nx,0:Np] = (latm.atm_data.x[1:Nx+1,0:Np] + latm.atm_data.x[0:Nx  ,1:Np+1] + latm.atm_data.x[1:Nx+1,1:Np+1])/3.0
    xc_tri4[0:Nx,0:Np] = (latm.atm_data.x[0:Nx  ,0:Np] + latm.atm_data.x[0:Nx  ,1:Np+1] + latm.atm_data.x[1:Nx+1,1:Np+1])/3.0

    pc_tri1[0:Nx,0:Np] = (latm.atm_data.p[0:Nx  ,0:Np] + latm.atm_data.p[1:Nx+1,0:Np  ] + latm.atm_data.p[0:Nx  ,1:Np+1])/3.0 
    pc_tri2[0:Nx,0:Np] = (latm.atm_data.p[0:Nx  ,0:Np] + latm.atm_data.p[1:Nx+1,0:Np  ] + latm.atm_data.p[1:Nx+1,1:Np+1])/3.0 
    pc_tri3[0:Nx,0:Np] = (latm.atm_data.p[1:Nx+1,0:Np] + latm.atm_data.p[0:Nx  ,1:Np+1] + latm.atm_data.p[1:Nx+1,1:Np+1])/3.0
    pc_tri4[0:Nx,0:Np] = (latm.atm_data.p[0:Nx  ,0:Np] + latm.atm_data.p[0:Nx  ,1:Np+1] + latm.atm_data.p[1:Nx+1,1:Np+1])/3.0
    
    
    # Slope of the two lines 1->3 and 2->4    
    A_1 = (pc_tri3[:,:]-pc_tri1[:,:])/(xc_tri3[:,:]-xc_tri1[:,:])
    A_2 = (pc_tri2[:,:]-pc_tri4[:,:])/(xc_tri2[:,:]-xc_tri4[:,:])
    
    # Intersection of the two lines
    latm.atm_data.xm[1:Nx+1,1:Np+1] = (A_1[:,:]*xc_tri1[:,:] - A_2[:,:]*xc_tri2[:,:] - pc_tri1[:,:] + pc_tri2[:,:])/(A_1[:,:] - A_2[:,:])
    latm.atm_data.pm[1:Nx+1,1:Np+1] = A_1[:,:]*(latm.atm_data.xm[1:Nx+1,1:Np+1] - xc_tri1[:,:]) + pc_tri1[:,:]
    
    
    # bottom,top boundary of pm
    for i in range(Nx):
        latm.atm_data.pm[i+1,0] = pA
        latm.atm_data.pm[i+1,Np+1]=(pB[i]+pB[i+1])/2.0
        
        
    # four corner points of pm
    #latm.atm_data.pm[0,0]=pA
    #latm.atm_data.pm[0,Np+1]=pB[0]            
    #latm.atm_data.pm[Nx+1,Np+1]=pB[Nx]
    #latm.atm_data.pm[Nx+1,0]=pA
        
    # left boundary of pm
    #latm.atm_data.pm[0,1:Np+1]=latm.atm_data.pm[Nx,1:Np+1]
    # right boundary of pm
    #latm.atm_data.pm[Nx+1,1:Np+1]=latm.atm_data.pm[1,1:Np+1]
    
    
    
    # bottom, top boundary of xm    
    latm.atm_data.xm[1:Nx+1,0] = latm.atm_data.xm[1:Nx+1,1] #(latm.atm_data.x[0:Nx,0] + latm.atm_data.x[1:Nx+1,0])/2.0
    latm.atm_data.xm[1:Nx+1,Np+1] = latm.atm_data.xm[1:Nx+1,Np] #(latm.atm_data.x[0:Nx,Np] + latm.atm_data.x[1:Nx+1,Np])/2.0
    
    # left right bonndary
    #for j in range(Np):
#        latm.atm_data.xm[0,j+1]=latm.atm_data.xm[Nx,j+1] - L
#        latm.atm_data.xm[Nx+1,j+1]=latm.atm_data.xm[1,j+1] + L
        
        
    #four corners of xm
    #latm.atm_data.xm[0,0]=latm.atm_data.xm[0,1]
    #latm.atm_data.xm[Nx+1,0]=latm.atm_data.xm[Nx+1,1]
    #latm.atm_data.xm[0,Np+1]=latm.atm_data.xm[0,Np]
    #latm.atm_data.xm[Nx+1,Np+1]=latm.atm_data.xm[Nx+1,Np]


    for j in range(Np+1):
        for i in range(Nx):
            latm.atm_data.normal[0,j*Nx+i]=-(latm.atm_data.p[i+1,j]-latm.atm_data.p[i,j])
            latm.atm_data.normal[1,j*Nx+i]=latm.atm_data.x[i+1,j]-latm.atm_data.x[i,j]
            norm=np.sqrt(latm.atm_data.normal[0,j*Nx+i]**2.0+latm.atm_data.normal[1,j*Nx+i]**2.0)
            latm.atm_data.hhori[j*Nx+i]=norm
            latm.atm_data.normal[0,j*Nx+i]=latm.atm_data.normal[0,j*Nx+i]/norm
            latm.atm_data.normal[1,j*Nx+i]=latm.atm_data.normal[1,j*Nx+i]/norm
            if j<Np:
                latm.atm_data.ac[i,j]=latm.atm_data.dx*np.abs(latm.atm_data.dp[i+1]+latm.atm_data.dp[i])/2.0
    
    # Find pm2 and pm3 meshes
    for j in range(Np):
        latm.atm_data.pm2[:,j+1]=(latm.atm_data.p[:,j]+latm.atm_data.p[:,j+1])/2.0
    latm.atm_data.pm2[:,0]=latm.atm_data.p[:,0]
    latm.atm_data.pm2[:,Np+1]=latm.atm_data.p[:,Np]
    
    for i in range(Nx):
        latm.atm_data.pm3[i+1,:]=(latm.atm_data.p[i+1,:]+latm.atm_data.p[i,:])/2.0
    latm.atm_data.pm3[0,:]=latm.atm_data.p[0,:]
    latm.atm_data.pm3[Nx+1,:]=latm.atm_data.p[Nx,:]

    
    
    
    
    #
    #  check this out!
    #
    latm.atm_data.pm[0,0:Np+2]=latm.atm_data.pm[Nx,0:Np+2]
    latm.atm_data.xm[0,0:Np+2]=latm.atm_data.xm[Nx,0:Np+2] - (xf-x0)
    latm.atm_data.pm[Nx+1,0:Np+2]=latm.atm_data.pm[1,0:Np+2]
    latm.atm_data.xm[Nx+1,0:Np+2]=latm.atm_data.xm[1,0:Np+2] + (xf-x0)
    
    latm.atm_data.pm3[0,0:Np+1]=latm.atm_data.pm3[Nx,0:Np+1]
    latm.atm_data.pm3[Nx+1,0:Np+1]=latm.atm_data.pm3[1,0:Np+1]
    
    latm.atm_data.xm3[0,0:Np+1]=latm.atm_data.xm3[Nx,0:Np+1] - (xf-x0)
    latm.atm_data.xm3[Nx+1,0:Np+1]=latm.atm_data.xm3[1,0:Np+1] + (xf-x0)
    
def random_qs(ep,Nx,Np):
	#print ep*np.random.normal(0.0,1.0 ,(Nx+2,Np+2))
	temp = np.zeros((Nx+2,Np+2))
	for i in range(1,13):
		temp = temp + 2.0*(np.random.uniform(0,1,(Nx+2,Np+2))-0.5)

	latm.atm_data.rqs = ep*temp/12.0
	
def random_bomb_place(bomb_num,ep_bomb,Nx,Np,bomb_flag):
	
	
	x_cell=np.zeros(bomb_num,'i',order='f')
	y_cell=np.zeros(bomb_num,'i',order='f')
	rand_val=np.zeros(bomb_num,'d',order='f')
	
	rx=2
	rp=1
	
	x_cell = np.random.randint(7-rx, size=bomb_num)
	y_cell = rp+np.random.randint(Np/4, size=bomb_num)
	rand_val=ep_bomb*np.random.normal(0, 1, bomb_num)

	#print latm.atm_data.xm[3*Nx/5,0]
	#print latm.atm_data.zm[3*Nx/5,Np], latm.atm_data.zm[Nx,Np]
	#print latm.atm_data.zm[3*Nx/5,3*Np/4], latm.atm_data.zm[Nx,3*Np/4]
	
	for i in range(bomb_num):
		latm.atm_data.a[2,3*Nx/5+x_cell[i]-rx:3*Nx/5+x_cell[i]+rx,Np-2-y_cell[i]-rp:Np-2-y_cell[i]+rp] = latm.atm_data.a[2,3*Nx/5+x_cell[i]-rx:3*Nx/5+x_cell[i]+rx,Np-2-y_cell[i]-rp:Np-2-y_cell[i]+rp] + rand_val[i]

	
def Average_Sol(Nx,Np):
	for i in range(Nx):
		latm.atm_data.a[:,i+1,:] = (latm.atm_data.a[:,i+1,:]+latm.atm_data.a[:,i,:])/2.0
		latm.atm_data.w[i+1,:] = (latm.atm_data.w[i+1,:]+latm.atm_data.w[i,:])/2.0
	
def CellK(Nx,Np,x,p,xm,pm):
    # Nx : Number of cells in x
    # Np : Number of cells in p
    # x : latm.data_data.x
    # p : latm.atm_data.p
    # xm : latm.atm_data.xm
    # pm : latm.atm_data.pm
    latm.atm_data.m = np.zeros((Nx,Np-1,4),'d', order='f')
    latm.atm_data.dux = np.zeros((Nx,Np-1),'d', order='f')
    
    latm.atm_data.ghost = np.zeros((2, Nx+2),'d',order = 'f')
    
    
    for i in range(Nx+2):
        x0 = latm.atm_data.xm[i,Np+1]
        p0 = latm.atm_data.pm[i,Np+1]
    
        if (i==0):
            x1 = latm.atm_data.xm2[0,Np+1] - latm.atm_data.dx
            p1 = latm.atm_data.pm2[Nx-1,Np+1]
            x2 = latm.atm_data.xm2[i,Np+1]
            p2 = latm.atm_data.pm2[i,Np+1]
        elif (i==Nx+1):            
            x1 = latm.atm_data.xm2[i-1,Np+1]
            p1 = latm.atm_data.pm2[i-1,Np+1]
            x2 = latm.atm_data.xm2[Nx,Np+1] + latm.atm_data.dx
            p2 = latm.atm_data.pm2[1,Np+1]
        else:
            x1 = latm.atm_data.xm2[i-1,Np+1]
            p1 = latm.atm_data.pm2[i-1,Np+1]
            x2 = latm.atm_data.xm2[i,Np+1]
            p2 = latm.atm_data.pm2[i,Np+1]
        
        m = (p2-p1)/(x2-x1)
        c = p1 - x1*m
        d = (x0+(p0-c)*m)/(1.0+m**2)
    
#        latm.atm_data.ghost[0,i] = 2.0*d - x0
#        latm.atm_data.ghost[1,i] = 2.0*d*m-p0+2.0*c
    
    
    
    
    for i in range(Nx):
        for j in range(Np-1):
            # M_{i,j+1/2} = [a b; c d]
            
            a=x[i+1,j+1]-x[i,j+1]
            b=p[i+1,j+1]-p[i,j+1]
            c=xm[i+1,j+2]-xm[i+1,j+1]
            d=pm[i+1,j+2]-pm[i+1,j+1]
            
            det=a*d-b*c
            latm.atm_data.m[i,j,0]=d/det
            latm.atm_data.m[i,j,1]=-b/det
            latm.atm_data.m[i,j,2]=-c/det
            latm.atm_data.m[i,j,3]=a/det
            
    latm.atm_data.coeff = np.zeros((Nx+1,Np-1,4),'d', order='f')
    latm.atm_data.cu = np.zeros((Nx+1,Np-1),'d', order='f')
    for i in range(Nx+1):
        for j in range(Np-1):
            latm.atm_data.coeff[i,j,0] = 0.25
            latm.atm_data.coeff[i,j,1] = (0.25*(4*pm[i,j+2]*x[i,j+1]-pm[i,j+2]*xm[i,j+1]-4*pm[i+1,j+2]*x[i,j+1]+pm[i+1,j+2]*xm[i,j+1]+4*xm[i+1,j+2]*p[i,j+1]-xm[i+1,j+2]*pm[i,j+1]-4*xm[i,j+2]*p[i,j+1]+xm[i,j+2]*pm[i,j+1]+3*xm[i,j+2]*pm[i+1,j+2]-3*xm[i+1,j+2]*pm[i,j+2])/(xm[i+1,j+2]*pm[i+1,j+1]-xm[i+1,j+2]*pm[i,j+2]-xm[i,j+2]*pm[i+1,j+1]+xm[i,j+2]*pm[i+1,j+2]+xm[i+1,j+1]*pm[i,j+2]-xm[i+1,j+1]*pm[i+1,j+2]))
            latm.atm_data.coeff[i,j,2] = (-0.25*(4*pm[i+1,j+1]*x[i,j+1]-pm[i+1,j+1]*xm[i,j+1]-4*pm[i+1,j+2]*x[i,j+1]+pm[i+1,j+2]*xm[i,j+1]+4*xm[i+1,j+2]*p[i,j+1]-xm[i+1,j+2]*pm[i,j+1]-4*xm[i+1,j+1]*p[i,j+1]+xm[i+1,j+1]*pm[i,j+1]+3*xm[i+1,j+1]*pm[i+1,j+2]-3*xm[i+1,j+2]*pm[i+1,j+1])/(xm[i+1,j+2]*pm[i+1,j+1]-xm[i+1,j+2]*pm[i,j+2]-xm[i,j+2]*pm[i+1,j+1]+xm[i,j+2]*pm[i+1,j+2]+xm[i+1,j+1]*pm[i,j+2]-xm[i+1,j+1]*pm[i+1,j+2]))
            latm.atm_data.coeff[i,j,3] = (0.25*(4*pm[i+1,j+1]*x[i,j+1]-pm[i+1,j+1]*xm[i,j+1]-4*pm[i,j+2]*x[i,j+1]+pm[i,j+2]*xm[i,j+1]+4*xm[i,j+2]*p[i,j+1]-xm[i,j+2]*pm[i,j+1]-4*xm[i+1,j+1]*p[i,j+1]+xm[i+1,j+1]*pm[i,j+1]+3*xm[i+1,j+1]*pm[i,j+2]-3*xm[i,j+2]*pm[i+1,j+1])/(xm[i+1,j+2]*pm[i+1,j+1]-xm[i+1,j+2]*pm[i,j+2]-xm[i,j+2]*pm[i+1,j+1]+xm[i,j+2]*pm[i+1,j+2]+xm[i+1,j+1]*pm[i,j+2]-xm[i+1,j+1]*pm[i+1,j+2]))
            
            
def Average_Sol(Nx,Np):
	for i in range(Nx):
		latm.atm_data.a[:,i+1,:] = (latm.atm_data.a[:,i+1,:]+latm.atm_data.a[:,i+2,:])/2.0
		latm.atm_data.w[i+1,:] = (latm.atm_data.w[i+1,:]+latm.atm_data.w[i+2,:])/2.0
		

###########################################################################
###########################################################################
def CellK2(Nx,Np,xm,pm):
    # Nx : Number of cells in x
    # Np : Number of cells in p
    # xm : latm.atm_data.xm
    # pm : latm.atm_data.pm
    latm.atm_data.m2 = np.zeros((Nx,Np,4),'d', order='f')    
    for i in range(Nx):
        for j in range(Np):
            # M_{i,j+1/2} = [a b; c d]
            

            a=xm[i+2,j+1]-xm[i,j+1]
            b=pm[i+2,j+1]-pm[i,j+1]
            c=xm[i+1,j+2]-xm[i+1,j]
            d=pm[i+1,j+2]-pm[i+1,j]
			
            
            det=a*d-b*c
            latm.atm_data.m2[i,j,0]=d/det
            latm.atm_data.m2[i,j,1]=-b/det
            latm.atm_data.m2[i,j,2]=-c/det
            latm.atm_data.m2[i,j,3]=a/det
###########################################################################    
###########################################################################    
def CellK3(Nx,Np,x,p,xm,pm,bc):
    # Nx : Number of cells in x
    # Np : Number of cells in p
    # x : latm.data_data.x
    # p : latm.atm_data.p
    # xm : latm.atm_data.xm
    # pm : latm.atm_data.pm
    # bc : 0 - periodic, 1 - Neumann.
    latm.atm_data.m3 = np.zeros((Nx+1,Np,4),'d', order='f')    
    for i in range(Nx+1):
        for j in range(Np):
            # M_{i,j+1/2} = [a b; c d]
            
            a=xm[i+1,j+1]-xm[i,j+1]
            b=pm[i+1,j+1]-pm[i,j+1]
            c=x[i,j+1]-x[i,j]
            d=p[i,j+1]-p[i,j]
            
            det=a*d-b*c
            latm.atm_data.m3[i,j,0]=d/det
            latm.atm_data.m3[i,j,1]=-b/det
            latm.atm_data.m3[i,j,2]=-c/det
            latm.atm_data.m3[i,j,3]=a/det
            
            
    latm.atm_data.coeff3 = np.zeros((Nx+1,Np+1,4),'d', order='f')
    
    for i in range(Nx+1):
        for j in range(Np+1):
            latm.atm_data.coeff3[i,j,0] = 0.25
            latm.atm_data.coeff3[i,j,1] =  -(xm[i,j]*pm[i+1,j+1] - pm[i,j]*xm[i+1,j+1] - xm[i,j]*pm[i,j+1] + pm[i,j]*xm[i,j+1] - 3*xm[i+1,j+1]*pm[i,j+1] + 3*pm[i+1,j+1]*xm[i,j+1] - 4*pm[i+1,j+1]*x[i,j] + 4*xm[i+1,j+1]*p[i,j] + 4*pm[i,j+1]*x[i,j] - 4*xm[i,j+1]*p[i,j])/(4*(xm[i+1,j]*pm[i+1,j+1] - pm[i+1,j]*xm[i+1,j+1] - xm[i+1,j]*pm[i,j+1] + pm[i+1,j]*xm[i,j+1] + xm[i+1,j+1]*pm[i,j+1] - pm[i+1,j+1]*xm[i,j+1]))
            latm.atm_data.coeff3[i,j,2] =  (xm[i,j]*pm[i+1,j] - pm[i,j]*xm[i+1,j] - xm[i,j]*pm[i,j+1] + pm[i,j]*xm[i,j+1] - 3*xm[i+1,j]*pm[i,j+1] + 3*pm[i+1,j]*xm[i,j+1] - 4*pm[i+1,j]*x[i,j] + 4*xm[i+1,j]*p[i,j] + 4*pm[i,j+1]*x[i,j] - 4*xm[i,j+1]*p[i,j])/(4*(xm[i+1,j]*pm[i+1,j+1] - pm[i+1,j]*xm[i+1,j+1] - xm[i+1,j]*pm[i,j+1] + pm[i+1,j]*xm[i,j+1] + xm[i+1,j+1]*pm[i,j+1] - pm[i+1,j+1]*xm[i,j+1]))
            latm.atm_data.coeff3[i,j,3] =  -(xm[i,j]*pm[i+1,j] - pm[i,j]*xm[i+1,j] - xm[i,j]*pm[i+1,j+1] + pm[i,j]*xm[i+1,j+1] - 3*xm[i+1,j]*pm[i+1,j+1] + 3*pm[i+1,j]*xm[i+1,j+1] - 4*pm[i+1,j]*x[i,j] + 4*xm[i+1,j]*p[i,j] + 4*pm[i+1,j+1]*x[i,j] - 4*xm[i+1,j+1]*p[i,j])/(4*(xm[i+1,j]*pm[i+1,j+1] - pm[i+1,j]*xm[i+1,j+1] - xm[i+1,j]*pm[i,j+1] + pm[i+1,j]*xm[i,j+1] + xm[i+1,j+1]*pm[i,j+1] - pm[i+1,j+1]*xm[i,j+1]))
 #-(a1*c2 - a2*c1 - a1*d2 + a2*d1 - 3*c1*d2 + 3*c2*d1 - 4*c2*x + 4*c1*y + 4*d2*x - 4*d1*y)/(4*(b1*c2 - b2*c1 - b1*d2 + b2*d1 + c1*d2 - c2*d1))
 # (a1*b2 - a2*b1 - a1*d2 + a2*d1 - 3*b1*d2 + 3*b2*d1 - 4*b2*x + 4*b1*y + 4*d2*x - 4*d1*y)/(4*(b1*c2 - b2*c1 - b1*d2 + b2*d1 + c1*d2 - c2*d1))
 #-(a1*b2 - a2*b1 - a1*c2 + a2*c1 - 3*b1*c2 + 3*b2*c1 - 4*b2*x + 4*b1*y + 4*c2*x - 4*c1*y)/(4*(b1*c2 - b2*c1 - b1*d2 + b2*d1 + c1*d2 - c2*d1))
#    latm.atm_data.coeff3[:,0,:]=0.25
#   latm.atm_data.coeff3[:,Np,:]=0.25
    
###########################################################################    
###########################################################################    

def coeff_r(Nx,Np,xm,xm3,pm,pm2):
    # Nx : Number of cells in x
    # Np : Number of cells in p
    # xm : latm.atm_data.xm
    # xm3 : latm.atm_data.xm3
    # pm : latm.atm_data.pm
    # pm2 : latm.atm_data.pm2
    
    latm.atm_data.rs=np.zeros((Nx,Np+1,2),'d', order='f')
    for i in range(Nx):
        for j in range(Np+1):
            x1=xm[i+1,j+1]
            x2=xm[i+1,j]
            x3=xm3[i+1,j]
            if( (x3 < max(x1,x2)) and (x3 >min(x1,x2))):
                latm.atm_data.rs[i,j,1]=np.absolute(x3-x1)
                latm.atm_data.rs[i,j,0]=np.absolute(x3-x2)
            else:
                latm.atm_data.rs[i,j,1]=np.absolute(x3-x1)
                latm.atm_data.rs[i,j,0]=np.absolute(x3-x2)
                if (latm.atm_data.rs[i,j,0] >= latm.atm_data.rs[i,j,1]):
                    latm.atm_data.rs[i,j,0]=-latm.atm_data.rs[i,j,0]
                else:
                    latm.atm_data.rs[i,j,1]=-latm.atm_data.rs[i,j,1]
            if (np.absolute(latm.atm_data.rs[i,j,0]+latm.atm_data.rs[i,j,1]) < 1E-8):
                latm.atm_data.rs[i,j,0]=1/2.0
                latm.atm_data.rs[i,j,1]=1/2.0
                
    latm.atm_data.rw=np.zeros((Nx+1,Np,2),'d', order='f')
    for i in range(Nx+1):
        for j in range(Np):
            p1=pm[i+1,j+1]
            p2=pm[i,j+1]
            p3=pm2[i,j+1]
            if ( (p3<max(p1,p2)) and (p3>min(p1,p2)) ):
                latm.atm_data.rw[i,j,1]=np.absolute(p3-p1)
                latm.atm_data.rw[i,j,0]=np.absolute(p3-p2)
            else:
                latm.atm_data.rw[i,j,1]=np.absolute(p3-p1)
                latm.atm_data.rw[i,j,0]=np.absolute(p3-p2)
                if (latm.atm_data.rw[i,j,0] >= latm.atm_data.rw[i,j,1]):
                    latm.atm_data.rw[i,j,0]=-latm.atm_data.rw[i,j,0]
                else:
                    latm.atm_data.rw[i,j,1]=-latm.atm_data.rw[i,j,1]
            if (np.absolute(latm.atm_data.rw[i,j,0]+latm.atm_data.rw[i,j,1]) < 1E-8):
             latm.atm_data.rw[i,j,0]=1/2.0
             latm.atm_data.rw[i,j,1]=1/2.0

## Compute $\frac{\partial }{\partial x}\int_{p_A}^{p_B} u dx$
def CompCond(Nx,Np,dp,dx,pm,u):
    ld=np.zeros((Nx),'d',order='f')
    for i in range(Nx):
        s1=0.0
        s2=0.0
        if i==Nx-1:
            for j in range(Np):
                s2=s2+u[i+2,j+1]*(dp[i+1]+dp[0])/2.0
                s1=s1+u[i,j+1]*(dp[i-1]+dp[i])/2.0
        elif i==0:
            for j in range(Np):
                s2=s2+u[i+2,j+1]*(dp[i+1]+dp[i+2])/2.0
                s1=s1+u[i,j+1]*(dp[i]+dp[i])/2.0
        else:
            for j in range(Np):
                s2=s2+u[i+2,j+1]*(dp[i+1]+dp[i+2])/2.0
                s1=s1+u[i,j+1]*(dp[i-1]+dp[i])/2.0

        ld[i]=(s2-s1)/(2*dx)/(pm[i+1,Np+1]-pm[i+1,0])
    return ld

## Compute $ \int_{p_A}^{p_B} u dx$
def CompCond2(Nx,Np,dp,dx,pm,u):
    pld=np.zeros((Nx),'d',order='f')
    for i in range(Nx):
        pld[i] = 0.0
        for j in range(Np):
            pld[i] = pld[i] + u[i+1,j+1]*(dp[i+1]+dp[i])/2.0        
        pld[i] = pld[i]
    return pld    
    
            
# latm.time_method.rk2source(method)
# Compute the Runge Kutta 2 method on all the mesh for the Shallow water system with coriolis force and source term.
# /!\ The fonction data_f90 has to re run before calling this function. /!\
#
# INPUT :
#    method = integer 0 or 1
#
# OUTPUT :
#    Q
def RK2Source():
    latm.time_method.rk2source()

# latm.time_method.rk4source(method)
# Compute the Runge Kutta 4 method on all the mesh for the Shallow water system with coriolis force and source term.
# /!\ The fonction data_f90 has to re run before calling this function. /!\
#
# INPUT :
#    method = integer 0 or 1
#
# OUTPUT :
#    Q
def RK4Source():
    latm.time_method.rk4source()
    
def Euler1():
    latm.time_method.euler1()
    
def Euler2():
    latm.time_method.euler2()

def Rk4_2():
    latm.time_method.rk4_2()


# Save data
# Save t,U,V in the file open as FILE_ID
def save_txt(DIR,t,Nx,Np,X,P,T,q,u,w):
    FILE_ID =  open(DIR+str(latm.atm_data.nt)+'.data', 'w')    
    FILE_ID.write(str(t) +'\t' + str(Nx)  + '\t' +str(Np) + '\n')
    for i in range(Nx+2):
        for j in range(Np+2):
            FILE_ID.write(str(X[i,j]) + '\t' + str(P[i,j]) + '\t' \
                + str(T[i,j]) + '\t' + str(q[i,j]) + '\t' \
                + str(u[i,j]) + '\t' + str(w[i,j]) + '\n')
    FILE_ID.close()
    
def print_mesh(Nx,Np,P):
    FILE_ID =  open('mesh_p.data', 'w')    
    FILE_ID.write(str(Nx)  + '\t' +str(Np) + '\n')
    for i in range(Nx+2):
        for j in range(Np+2):
            FILE_ID.write(str(P[i,j]) + '\n')
    FILE_ID.close()
    
def fz(p,deltaT, g, p0, T0):
    R = 287.0
    C = R*(T0-deltaT)*np.log(p0) + R*deltaT
    return (-R*(T0-deltaT)*np.log(p) - R*deltaT/p0*p +C)/g


# read text file
# DIR = file path
def read_txt(FileName):
    
    f1 = file(FileName, "r")
    temp_ele = [line.split() for line in f1]
    f1.close()
    start_line = map(float,temp_ele[0])
    t = start_line[0]
    Nx = int(start_line[1])
    Np = int(start_line[2])
        
    #initialize
    X = np.zeros((Nx+2,Np+2),'d',order='f')
    P = np.zeros((Nx+2,Np+2),'d',order='f')
    T = np.zeros((Nx+2,Np+2),'d',order='f')
    q = np.zeros((Nx+2,Np+2),'d',order='f')
    u = np.zeros((Nx+2,Np+2),'d',order='f')
    w = np.zeros((Nx+2,Np+2),'d',order='f')    
    
    l=1
    for i in range(Nx+2):
        for j in range(Np+2):
            temp = map(float,temp_ele[l])
            X[i,j] = temp[0]
            P[i,j] = temp[1]
            T[i,j] = temp[2]
            q[i,j] = temp[3]
            u[i,j] = temp[4]
            w[i,j] = temp[5]
            l=l+1
    return [t,Nx,Np,X,P,T,q,u,w]
	
