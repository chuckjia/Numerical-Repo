# PLoting libary

import matplotlib.pyplot as plt
import numpy as np

def plot_periodic(fig,x0,xf,Nx,Np,T,q,t,nt,xm,zm,colors1,colors2,title1,title2,name,SAVEFIG,PLOTFIG):
	Tper=np.zeros((3*Nx,Np+2),'d',order='f')
	qper=np.zeros((3*Nx,Np+2),'d',order='f')
	xmper=np.zeros((3*Nx,Np+2),'d',order='f')
	zmper=np.zeros((3*Nx,Np+2),'d',order='f')
	
	Tper[0:Nx,:]=T[1:Nx+1,:]
	Tper[Nx:2*Nx,:]=T[1:Nx+1,:]
	Tper[2*Nx:3*Nx,:]=T[1:Nx+1,:]
	
	qper[0:Nx,:]=q[1:Nx+1,:]
	qper[Nx:2*Nx,:]=q[1:Nx+1,:]
	qper[2*Nx:3*Nx,:]=q[1:Nx+1,:]

	xmper[0:Nx,:]=xm[1:Nx+1,:]-(xf-x0)
	xmper[Nx:2*Nx,:]=xm[1:Nx+1,:]
	xmper[2*Nx:3*Nx,:]=xm[1:Nx+1,:]+(xf-x0)

	zmper[0:Nx,:]=zm[1:Nx+1,:]
	zmper[Nx:2*Nx,:]=zm[1:Nx+1,:]
	zmper[2*Nx:3*Nx,:]=zm[1:Nx+1,:]
	
	#print "toto",xm,zm,xmper,zmper

	plot_contour(fig,Tper,qper,t,nt,xmper,zmper,xmper,zmper,colors1,colors2,title1,title2,name,SAVEFIG,PLOTFIG)
			
#Figures With overlappings.
	     
def plot_2d_single_over(fig,val1,val2,t,nt,x1,x2,title1,title2,name,SAVEFIG,PLOTFIG):
	#plt.figure(fig,figsize=(40,40))
	plt.figure(fig)	
	plt.clf()
	
	plt.plot(x1,val1,'k--', label='before the projection')
	plt.plot(x2,val2,'k', label='after the projection')
	plt.legend(loc='upper right')
	#plt.axis(axis)
	NewTitle = title1 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	#plt.xlabel('x')


	if (SAVEFIG==1):
		if (nt<10):
			filename = name + '_000000000' + str(nt) + '.png'
		elif (nt<100):
			filename = name + '_00000000' + str(nt) + '.png'
		elif (nt<1000):
			filename = name + '_0000000' + str(nt) + '.png'
		elif (nt<10000):
			filename = name + '_000000' + str(nt) + '.png'
		elif (nt<100000):
			filename = name + '_00000' + str(nt) + '.png'
		elif (nt<1000000):
			filename = name + '_0000' + str(nt) + '.png'
		elif (nt<10000000):
			filename = name + '_000' + str(nt) + '.png'
		elif (nt<100000000):
			filename = name + '_00' + str(nt) + '.png'
		elif (nt<1000000000):
			filename = name + '_0' + str(nt) + '.png'
		else:
			filename = name + '_' + str(nt) + '.png'
		plt.savefig(filename)

	if (PLOTFIG==1):
		#plt.ion()
		plt.draw()

			

# plot_error(fig,errh,t)
def plot_error(fig,err1,err2,err3,tv,title,leg1,leg2,leg3):
	plt.figure(fig)
	plt.clf()
	plt.plot(tv,err1,'b',tv,err2,'r',tv,err3,'g')
	plt.title(title)
	plt.legend( (leg1, leg2, leg3) )
	filename='IMG/'+title	
	plt.savefig(filename)
	#plt.show()

# plot_contour(fig,h,t,nt,X,Y,title,name,SAVEFIG,PLOTFIG)
# plot or save the plotting (in contour) of the variable h on the mesh [X,Y]
#
# INPUT :
#	fig = integer > 0 which is the number for the figure
#	h = double of dimension (Nx,Ny)
#	t = double with is the time at what the variable h is computed
#	nt = integer which the number of the plot
#	X = double of dimension Nx
#	Y = double of dimension Ny
#	title = chain of character title of the plot
#	name = chain of character name of the file in what the figure will be saved
#	SAVEFIG = integer, 1 figure saved, 0 figure not saved
#	PLOTFIG = integer, ! figure ploted, 0 figure not ploted
def plot_contour(fig,T,q,t,nt,X1,P1,X2,P2,colors1,colors2,title1,title2,name,SAVEFIG,PLOTFIG):
	#plt.figure(fig,figsize=(10,5))
	plt.figure(fig)
	plt.clf()
	
	plt.subplot(211)
	#plt.axis("equal")
	if len(colors1)==0:
		plt.contourf(X1,P1,T)
	else:
		plt.contourf(X1,P1,T,colors1)
	plt.colorbar(extend="both",format="%.2e")
	#plt.ylim(pB,p0)
	NewTitle = title1 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	
	plt.subplot(212)
	#plt.axis("equal")
	if len(colors2)==0:
		plt.contourf(X2,P2,q)
	else:
		plt.contourf(X2,P2,q,colors2)
	plt.colorbar(extend="both",format="%.2e")
	#plt.ylim(pB,p0)
	NewTitle = title2 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	
	
	# To set y axis
	# plt.ylim(-2,2)
	# See http://matplotlib.org/api/pyplot_api.html
	
	
	
	#plt.colorbar.set_clim(vmin=0,vmax=2.0) 
	if (SAVEFIG==1):
		if (nt<10):
			filename = name + '_000000000' + str(nt) + '.png'
		elif (nt<100):
			filename = name + '_00000000' + str(nt) + '.png'
		elif (nt<1000):
			filename = name + '_0000000' + str(nt) + '.png'
		elif (nt<10000):
			filename = name + '_000000' + str(nt) + '.png'
		elif (nt<100000):
			filename = name + '_00000' + str(nt) + '.png'
		elif (nt<1000000):
			filename = name + '_0000' + str(nt) + '.png'
		elif (nt<10000000):
			filename = name + '_000' + str(nt) + '.png'
		elif (nt<100000000):
			filename = name + '_00' + str(nt) + '.png'
		elif (nt<1000000000):
			filename = name + '_0' + str(nt) + '.png'
		else:
			filename = name + '_' + str(nt) + '.png'
		plt.savefig(filename)
	if (PLOTFIG!=0):
		plt.ion()
		plt.draw()




def plot_contour_single(fig,T,t,nt,X1,P1,colors1,title1,name,SAVEFIG,PLOTFIG):
	#plt.figure(fig,figsize=(10,5))
	plt.figure(fig)
	plt.clf()
	
	#plt.subplot(211)
	#plt.axis("equal")
	if len(colors1)==0:
		plt.contourf(X1,P1,T)
	else:
		plt.contourf(X1,P1,T,colors1)
	plt.colorbar(extend="both",format="%.2e")
	#plt.ylim(pB,p0)
	NewTitle = title1 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	
	#plt.subplot(212)
	##plt.axis("equal")
	#if len(colors2)==0:
		#plt.contourf(X2,P2,q)
	#else:
		#plt.contourf(X2,P2,q,colors2)
	#plt.colorbar(extend="both",format="%.2e")
	##plt.ylim(pB,p0)
	#NewTitle = title2 + " at t=" + str(t)
	#plt.title(NewTitle,fontsize=12)
	
	
	# To set y axis
	# plt.ylim(-2,2)
	# See http://matplotlib.org/api/pyplot_api.html
	
	
	
	#plt.colorbar.set_clim(vmin=0,vmax=2.0) 
	if (SAVEFIG==1):
		if (nt<10):
			filename = name + '_000000000' + str(nt) + '.png'
		elif (nt<100):
			filename = name + '_00000000' + str(nt) + '.png'
		elif (nt<1000):
			filename = name + '_0000000' + str(nt) + '.png'
		elif (nt<10000):
			filename = name + '_000000' + str(nt) + '.png'
		elif (nt<100000):
			filename = name + '_00000' + str(nt) + '.png'
		elif (nt<1000000):
			filename = name + '_0000' + str(nt) + '.png'
		elif (nt<10000000):
			filename = name + '_000' + str(nt) + '.png'
		elif (nt<100000000):
			filename = name + '_00' + str(nt) + '.png'
		elif (nt<1000000000):
			filename = name + '_0' + str(nt) + '.png'
		else:
			filename = name + '_' + str(nt) + '.png'
		plt.savefig(filename)
	if (PLOTFIG!=0):
		plt.ion()
		plt.draw()





def plot_curve(fig,T,q,t,nt,X1,P1,X2,P2,colors1,colors2,title1,title2,name,SAVEFIG,PLOTFIG):
	#plt.figure(fig,figsize=(10,5))
	plt.figure(fig)
	plt.clf()
	
	plt.subplot(211)
	#plt.axis("equal")
	if len(colors1)==0:
		plt.contour(X1,P1,T)
	else:
		plt.contour(X1,P1,T,colors1)
	plt.colorbar(extend="both",format="%.2e")
	#plt.ylim(pB,p0)
	NewTitle = title1 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	
	plt.subplot(212)
	#plt.axis("equal")
	if len(colors2)==0:
		plt.contour(X2,P2,q)
	else:
		plt.contour(X2,P2,q,colors2)
	plt.colorbar(extend="both",format="%.2e")
	#plt.ylim(pB,p0)
	NewTitle = title2 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	
	
	# To set y axis
	# plt.ylim(-2,2)
	# See http://matplotlib.org/api/pyplot_api.html
	
	
	
	#plt.colorbar.set_clim(vmin=0,vmax=2.0) 
	if (SAVEFIG==1):
		if (nt<10):
			filename = name + '_000000000' + str(nt) + '.png'
		elif (nt<100):
			filename = name + '_00000000' + str(nt) + '.png'
		elif (nt<1000):
			filename = name + '_0000000' + str(nt) + '.png'
		elif (nt<10000):
			filename = name + '_000000' + str(nt) + '.png'
		elif (nt<100000):
			filename = name + '_00000' + str(nt) + '.png'
		elif (nt<1000000):
			filename = name + '_0000' + str(nt) + '.png'
		elif (nt<10000000):
			filename = name + '_000' + str(nt) + '.png'
		elif (nt<100000000):
			filename = name + '_00' + str(nt) + '.png'
		elif (nt<1000000000):
			filename = name + '_0' + str(nt) + '.png'
		else:
			filename = name + '_' + str(nt) + '.png'
		plt.savefig(filename)
	if (PLOTFIG!=0):
		plt.ion()
		plt.draw()

def plot_quiver(fig,U,W,t,nt,title,name,SAVEFIG,PLOTFIG):
	#plt.figure(fig,figsize=(10,5))
	plt.figure(fig)
	plt.clf()
	
	plt.quiver(U,W)
	NewTitle = title + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
		
	
	#plt.colorbar.set_clim(vmin=0,vmax=2.0) 
	if (SAVEFIG==1):
		if (nt<10):
			filename = name + '_000000000' + str(nt) + '.png'
		elif (nt<100):
			filename = name + '_00000000' + str(nt) + '.png'
		elif (nt<1000):
			filename = name + '_0000000' + str(nt) + '.png'
		elif (nt<10000):
			filename = name + '_000000' + str(nt) + '.png'
		elif (nt<100000):
			filename = name + '_00000' + str(nt) + '.png'
		elif (nt<1000000):
			filename = name + '_0000' + str(nt) + '.png'
		elif (nt<10000000):
			filename = name + '_000' + str(nt) + '.png'
		elif (nt<100000000):
			filename = name + '_00' + str(nt) + '.png'
		elif (nt<1000000000):
			filename = name + '_0' + str(nt) + '.png'
		else:
			filename = name + '_' + str(nt) + '.png'
		plt.savefig(filename)
	if (PLOTFIG!=0):
		plt.ion()
		plt.draw()


def plot_2d(fig,val1,val2,t,nt,x1,x2,title1,title2,name,SAVEFIG,PLOTFIG):
	#plt.figure(fig,figsize=(40,40))
	plt.figure(fig)	
	plt.clf()
	
	plt.subplot(211)
	plt.plot(x1,val1)
	#plt.axis(axis)
	NewTitle = title1 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	#plt.xlabel('x')
	
	plt.subplot(212)
	plt.plot(x2,val2)
	#plt.axis(axis)
	NewTitle = title2 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	#plt.xlabel('x')

	if (SAVEFIG==1):
		if (nt<10):
			filename = name + '_000000000' + str(nt) + '.png'
		elif (nt<100):
			filename = name + '_00000000' + str(nt) + '.png'
		elif (nt<1000):
			filename = name + '_0000000' + str(nt) + '.png'
		elif (nt<10000):
			filename = name + '_000000' + str(nt) + '.png'
		elif (nt<100000):
			filename = name + '_00000' + str(nt) + '.png'
		elif (nt<1000000):
			filename = name + '_0000' + str(nt) + '.png'
		elif (nt<10000000):
			filename = name + '_000' + str(nt) + '.png'
		elif (nt<100000000):
			filename = name + '_00' + str(nt) + '.png'
		elif (nt<1000000000):
			filename = name + '_0' + str(nt) + '.png'
		else:
			filename = name + '_' + str(nt) + '.png'
		plt.savefig(filename)

	if (PLOTFIG==1):
		#plt.ion()
		plt.draw()

def plot_2d_single(fig,val1,t,nt,x1,title1,name,SAVEFIG,PLOTFIG):
	#plt.figure(fig,figsize=(40,40))
	plt.figure(fig)	
	plt.clf()
	
	plt.plot(x1,val1)
	#plt.axis(axis)
	NewTitle = title1 + " at t=" + str(t)
	plt.title(NewTitle,fontsize=12)
	#plt.xlabel('x')

	if (SAVEFIG==1):
		if (nt<10):
			filename = name + '_000000000' + str(nt) + '.png'
		elif (nt<100):
			filename = name + '_00000000' + str(nt) + '.png'
		elif (nt<1000):
			filename = name + '_0000000' + str(nt) + '.png'
		elif (nt<10000):
			filename = name + '_000000' + str(nt) + '.png'
		elif (nt<100000):
			filename = name + '_00000' + str(nt) + '.png'
		elif (nt<1000000):
			filename = name + '_0000' + str(nt) + '.png'
		elif (nt<10000000):
			filename = name + '_000' + str(nt) + '.png'
		elif (nt<100000000):
			filename = name + '_00' + str(nt) + '.png'
		elif (nt<1000000000):
			filename = name + '_0' + str(nt) + '.png'
		else:
			filename = name + '_' + str(nt) + '.png'
		plt.savefig(filename)

	if (PLOTFIG==1):
		#plt.ion()
		plt.draw()


# 3D plot with Gnuplot
def plot_3D(h,xm,ym):
	g = gplt.Gnuplot(debug=1)
	g.splot(gplt.GridData(h,xm,ym, binary=0))
	raw_input('Please press return to continue...\n')
	
