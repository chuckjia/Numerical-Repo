import numpy as np

# Resolution  with source term of Shallow Water equations with Coriolis force :
# dQ/dt + dF(Q)/dx + dG(Q)/dy + C(Q) = S(t)
#
# With :
# 	Q = transpose(h,U,V) (data)
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

# Data
CLEAN = 'yes'
SOURCE = 0
DEBUG = 0
SAVEFIG = 1
FIG_NAME = 'atmoshpere'
PLOTFIG = 0 # 1: plot, 0 : no plot
PLOT_TYPE = 2 # 1: plot_all, 2: plot_contour, 3 : plot 3D mayavi, 4: plot 2D
titleT = "T"
titleq = "q"
#colors1=np.linspace(250.0,320.0,10)
#colors2=np.linspace(0.0,0.00026,10)
#colors2=np.linspace(0.0,0.028,10)
#colors2=np.linspace(0.07,0.14,10)
colors5=np.linspace(255.0,308.0,30)
colors7=np.linspace(255.0,330.0,30)
colors8=np.linspace(8.0*10**(-3),2.4*10**(-2),12)
#colors2=np.linspace(0.0,0.00026,10)
#colors2=np.linspace(0.0,0.028,10)
colors2=np.linspace(0.0015,0.018,30)
colors1=[]
colors2b=[]
colors3=[]
colors4=[]
#method = 4 # 1 = center , 2 = Centra-Upwind
time_method = 5 # 2=RK2, 4=RK4
scheme= 1
th1=1.0
th2=2.0
SAVE2TXT_DATA= 0
PERIODE_FIG = 1000
PERIODE_DATA = 10*1E10
PLOT_ERROR=0
PLOT_energy=1
PERIODE_ERROR=0
PHYSIC = 1
# bc=1 : Neumann-Neumann
# bc=2 : Dirichlet-Neumann
bc=2
# Data for the Mesh
p0=1000.0
#pA=10.0
pA=250.0
x0=0.0
xf=75000.0
Nx=100#2500
Np=100#45
u0=0.02
g=9.8
deltaT=50.0
T0=300.0
#ep = 0.003
ep = 0.0
# Time discretization
nt0=0
t0=0.0
dt=0.5
tf=t0+60000#59200.0
pi=np.pi
AVG1 = 18
energy_count = 1
Ep_dump=0.0#1.5*1E-4


bombing=0
bomb_num=48
ep_bomb=5.0*np.sqrt(dt)
bomb_period=1
bomb_lenght=1
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
				pbt[i,j] = 1000.0 -pA*exp(-((x[i,j]-(xf/2.0))**2)/6000.0**2)				
	else:
		for i in range(len(x)):
			pbt[i] = 1000.0 -pA*exp(-((x[i]-(xf/2.0))**2)/6000.0**2)						 
	return pbt

# Exact Solution
def T_init(x,p):
	return T0-(1.0-p/p0)*deltaT #+ cos(p0/pB_ex(x)*pi)*10.0
	#pbt=np.array(x)
	#for i in range(len(x[:,0])):
		#for j in range(len(x[0,:])):
			#if (cos(pi*x[i,j]/xf)/100.0 <=  0):
				#pbt[i,j]= T0-(1.0-p[i,j]/p0)*deltaT - 10.0*cos(p0*2*pi)
			#else:
				#pbt[i,j]= T0-(1.0-p[i,j]/p0)*deltaT
	#return pbt	
	
	
def T_bar(x,p):
	return T0-(1.0-p/p0)*deltaT
	
def q_init(x,p):
	pbt=np.array(x)
	for i in range(len(x[:,0])):
		for j in range(len(x[0,:])):
			if (x[i,j]/xf <=  0.2):
				pbt[i,j]= ((0.622*6.112*np.exp( (17.65*(T_bar(x[i,j],p[i,j])-273.15))/(T_bar(x[i,j],p[i,j])-29.65) ))/p[i,j])  - 0.0052 #+ (cos(2.5*pi*x[i,j]/xf)/100.0)*(p[i,j]/1000.0)
			else:
				pbt[i,j] = 0.0
				pbt[i,j]= ((0.622*6.112*np.exp( (17.65*(T_bar(x[i,j],p[i,j])-273.15))/(T_bar(x[i,j],p[i,j])-29.65) ))/p[i,j])  - 0.0052
			if (pbt[i,j] <= 0):
				pbt[i,j] = 0.0
	return pbt	

	#for i in range(len(x[:,0])):
		#for j in range(len(x[0,:])):
			#if (cos(pi*x[i,j]/xf)/100.0 <=  0):
				#pbt[i,j]= 0.6220/p[i,j]*2.53*10.0**(8.0)*exp( -5.43*10**3/T_bar(x[i,j],p[i,j] ) ) - 0.00001 + cos(pi*x[i,j]/xf)/1000.0
			#else:
				#pbt[i,j]= 0.6220/p[i,j]*2.53*10.0**(8.0)*exp( -5.43*10**3/T_bar(x[i,j],p[i,j] ) ) - 0.00001
			#if (pbt[i,j] <= 0):
				#pbt[i,j] = 0.0
	#return pbt		
	
def q_s(x,p):
	return ((0.622*6.112*np.exp( (17.65*(T_init(x,p)-273.15))/(T_init(x,p)-29.65) ))/p)
	#return 0.6220/p*2.53*10.0**(8.0)*exp( -5.43*10**3/T_init(x,p) )
	
def u_init(x,p):
	return 7.5+cos(p*pi/p0)*cos(2*pi*x/xf)

def wb(x,t):
	return w_ex(x,pB_ex(x),t)

def w_init(x,p):
	return 2.0*sin(2*pi*x/xf)*2.0*pi/xf*(sin(p*pi/p0)-sin(pA*pi/p0))*p0/pi#10.0**(-10)*x#

def Sw(x,p,t):
	return 0.0

def ST(x,p,t):
	return 0.0

def Sq(x,p,t):
	return 0.0

def Su(x,p,t):
	return 0.0
