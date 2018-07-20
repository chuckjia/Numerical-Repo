#!/usr/bin/env python2.7

import numpy
from mlplot import *
import matplotlib.pyplot as plt
import struct
import scipy.io as sio

DIR='run_physic4e3/'
tf=600000
dt=0.5
Nx=100

convert_text=1

Nit=int(tf/1/dt)

L2_T=np.zeros((Nx,Nit) ,'d')
#L2_q=np.zeros((Nx,Nit) ,'d')
L2_u=np.zeros((Nx,Nit) ,'d')
#L2_w=np.zeros((Nx,Nit) ,'d')

#t_energy = np.zeros(Nit ,'d')

doublesize = struct.calcsize('d')

fbin_T=open(DIR+'T_profile.dat', 'rb')
#fbin_q=open(DIR+'q_profile.dat', 'rb')
fbin_u=open(DIR+'u_profile.dat', 'rb')
#fbin_w=open(DIR+'w_profile.dat', 'rb')


#fbin_time=open(DIR+'L2_time.dat', 'rb')

	

for i in range(Nit):
	
	for j in range(Nx):
		bin_T = fbin_T.read(doublesize)
		data_T=struct.unpack('d', bin_T)
		L2_T[j,i]=data_T[0]
		
		#bin_q = fbin_q.read(doublesize)
		#data_q=struct.unpack('d', bin_q)
		#L2_q[j,i]=data_q[0]
		
		bin_u = fbin_u.read(doublesize)
		data_u=struct.unpack('d', bin_u)
		L2_u[j,i]=data_u[0]
		
		#bin_w = fbin_w.read(doublesize)
		#data_w=struct.unpack('d', bin_w)
		#L2_w[j,i]=data_w[0]
	
	#bin_t = fbin_time.read(doublesize)
	#data_t=struct.unpack('d', bin_t)
	#t_energy[i]=data_t[0]


if convert_text==1:
	sio.savemat(DIR+'T_profile.mat', {'T_profile':L2_T})
	#sio.savemat(DIR+'q_profile.mat', {'q_profile':L2_q})
	sio.savemat(DIR+'u_profile.mat', {'u_profile':L2_u})
	#sio.savemat(DIR+'w_profile.mat', {'w_profile':L2_w})
	#sio.savemat(DIR+'time.mat', {'time':t_energy})


fbin_T.close()
#fbin_q.close()
fbin_u.close()
#fbin_w.close()
#fbin_time.close()



#plt.figure(1)	
#plt.clf()
#plt.plot(range(Nx),L2_T[:,0])
#plt.savefig(DIR+'TTTT.png')




