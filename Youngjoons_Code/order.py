#!/usr/bin/env python2.6

import matplotlib.pyplot as plt
import numpy as np

print "For Nx=50, Np=50, dt=1E-4, t0=10"
print 't= 10.001  log(dx)= -0.698970004336  log(dp)= -0.690196080029  log(errT)= -7.36960809314  log(errq)= -7.36960805342  log(erru)= -6.39956417854  log(errw)= -6.52150619537'

print "For Nx=100, Np=100, dt=1E-4, t0=10"
print 't= 10.001  log(dx)= -1.0  log(dp)= -0.995635194598  log(errT)= -7.82259494023  log(errq)= -7.82259488396  log(erru)= -6.55204794018  log(errw)= -6.81321357319'

print "For Nx=200, Np=200, dt=1E-4, t0=10"
print 't= 10.001  log(dx)= -1.30102999566  log(dp)= -1.29885307641  log(errT)= -8.2782517537  log(errq)= -8.27825169846  log(erru)= -6.70360047682  log(errw)= -7.00757086527'

logdx=np.zeros(3,'d')
logdx=[-0.698970004336,-1.0,-1.30102999566]
logerrT=[-7.36960809314,-7.82259494023,-8.2782517537]
logerrq=[-7.36960805342,-7.82259488396,-8.27825169846 ]
logerru=[-6.39956417854,-6.55204794018,-6.70360047682]
logerrw=[-6.52150619537,-6.81321357319,-7.00757086527]

print "Slopes"
print "for T", (logerrT[2]-logerrT[1])/(logdx[2]-logdx[1])
print "for q", (logerrq[2]-logerrq[1])/(logdx[2]-logdx[1])
print "for u", (logerru[2]-logerru[1])/(logdx[2]-logdx[1])
print "for w", (logerrw[2]-logerrw[1])/(logdx[2]-logdx[1])

plt.figure(1)	
plt.clf()
plt.plot(logdx,logerrT,'b-',logdx,logerrq,'r-',logdx,logerru,'g-',logdx,logerrw,'k-')
plt.legend( ('log10(errT)', 'log10(errq)', 'log10(erru)','log19(errw)') )
plt.xlabel('log10(dx)')
plt.title('Error in $L^2$ at t=10.001')
filename="Error_in_space"	
plt.savefig(filename)
