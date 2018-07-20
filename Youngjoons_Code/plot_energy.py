#!/usr/bin/env python2.7

import numpy
from mlplot import *
import matplotlib.pyplot as plt
import struct
import scipy.io as sio

DIR='run_physic4e3/'
tf=600000
dt=0.5

convert_text=1

Nit=int(tf/1/dt)-1

L2_T=np.zeros(Nit ,'d')
#L2_q=np.zeros(Nit ,'d')
L2_u=np.zeros(Nit ,'d')
#L2_w=np.zeros(Nit ,'d')
t_energy = np.zeros(Nit ,'d')

doublesize = struct.calcsize('d')

fbin_T=open(DIR+'L2_T.dat', 'rb')
#fbin_q=open(DIR+'L2_q.dat', 'rb')
fbin_u=open(DIR+'L2_u.dat', 'rb')
#fbin_w=open(DIR+'L2_w.dat', 'rb')
fbin_t=open(DIR+'L2_time.dat', 'rb')

if convert_text==1:
	f=open(DIR+'L2_energy.txt', 'w')
	

for i in range(Nit):
	
	bin_T = fbin_T.read(doublesize)
	data_T=struct.unpack('d', bin_T)
	L2_T[i]=data_T[0]
	
	#bin_q = fbin_q.read(doublesize)
	#data_q=struct.unpack('d', bin_q)
	#L2_q[i]=data_q[0]
	
	bin_u = fbin_u.read(doublesize)
	data_u=struct.unpack('d', bin_u)
	L2_u[i]=data_u[0]
	
	#bin_w = fbin_w.read(doublesize)
	#data_w=struct.unpack('d', bin_w)
	#L2_w[i]=data_w[0]
	
	bin_t = fbin_t.read(doublesize)
	data_t=struct.unpack('d', bin_t)
	t_energy[i]=data_t[0]
	

	
	#if convert_text==1:
		#T='%.4E' % L2_T[i]
		#q='%.4E' % L2_q[i]
		#u='%.4E' % L2_u[i]
		#w='%.4E' % L2_w[i]
		#t=t_energy[i]
		#f.write('||T||= '+T+'\t'+'||q||= '+q+'\t'+'||u||= '+u+'\t'+'||w||= '+w+'\t'+'t= '+str(t)+'\n')

if convert_text==1:
	sio.savemat(DIR+'L2_T.mat', {'T':L2_T})
	#sio.savemat(DIR+'L2_q.mat', {'q':L2_q})
	sio.savemat(DIR+'L2_u.mat', {'u':L2_u})
	#sio.savemat(DIR+'L2_w.mat', {'w':L2_w})
	sio.savemat(DIR+'L2_time.mat', {'time':t_energy})

fbin_T.close()
#fbin_q.close()
fbin_u.close()
#fbin_w.close()
fbin_t.close()

if convert_text==1:
	f.close()

#plt.figure(1)	
#plt.clf()

#plt.subplot(211)
#plt.plot(t_energy,L2_T)
##plt.plot(t_energy[3000:60000:10],L2_T[3000:60000:10])
##plt.xlabel('t')
#plt.title('L2 of T over time',fontsize=12)

#plt.subplot(212)
#plt.plot(t_energy,L2_q)
##plt.plot(t_energy[3000:60000:10],L2_q[3000:60000:10])
#plt.xlabel('t')
#plt.title('L2 of q over time',fontsize=12)

#plt.draw()
#plt.savefig(DIR+'L2_energy_T_q.png')	

#plt.figure(1)	
#plt.clf()

#plt.subplot(211)
#plt.plot(t_energy,L2_u)
##plt.plot(t_energy[3000:60000:10],L2_u[3000:60000:10])
#plt.xlabel('t')
#plt.title('L2 of u over time',fontsize=12)

#plt.subplot(212)
#plt.plot(t_energy,L2_w)
##plt.plot(t_energy[3000:60000:10],L2_w[3000:60000:10])
#plt.xlabel('t')
#plt.title('L2 of w over time',fontsize=12)

#plt.draw()
#plt.savefig(DIR+'L2_energy_u_w.png')


