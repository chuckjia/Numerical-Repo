import numpy as np

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== # 
# Meshing parameters
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== # 

Nx = 10
Np = 10



# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== # 
# Geometrical parameters
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== # 

# x direction

x0 = 0.
xL = 50000.

Dx = (xL - x0) / Nx

# p direction

pA = 200.

def pB(x):
    return 1000. - 200. * np.exp(- (x - 25000.) ** 2. / 9000000.)

def Dp(x):
    return (pB(x) - pA) / Np




