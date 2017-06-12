from Constants import *
import numpy as np
from Carbon.QuickTime import suppressDither

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Initial conditions: wrapper functions
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

# Initial condition for T

def TInitial(x, p):
    return TManufactured2(x, p, 0)

# ----- ----- ----- ----- ----- ----- #

# Initial condition for q

def qInitial(x, p):
    return 0


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Source term: wrapper function
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

def calcS(T, q, x, p, t):
    """Calculates the value of the source term. 
       Return type: list of 2 elements."""
    return SourceManufactured2(T, q, x, p, t)


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Tests
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

# ----- Test 1 ----- #

def TInitialTest1(x, p, t):
    """The initial condition in the paper for T"""
    return 300 + (1 - p / p0) * 50

def qInitialTest1(x, p):
    """The initial condition in the paper for q"""
    T = TInitialTest1(x, p)
    qs = 0.622 * 6.112 * np.exp(17.67 * (T - 273.15) / (T - 29.65)) / p
    return qs - 0.0052


# ----- Test 2 ----- #

suppLeft = x0 + 4. * (xL - x0) / 9.
suppRight = x0 + 5. * (xL - x0) / 9.
pSec = (pB(0) - pA) / 9.
suppBottom = pA + 4. * pSec
suppTop = suppBottom + pSec

def TManufactured1(x, p, t):
    if x >= suppLeft and x <= suppRight and p >= suppBottom and p <= suppTop:
        xTemp = (2. * x - suppLeft - suppRight) / (suppRight - suppLeft)
        pTemp = (2. * p - suppBottom - suppTop) / (suppTop - suppBottom)
        return np.cos(t) \
             * np.exp(-1. / (1. - xTemp ** 2.)) \
             * np.exp(-1. / (1. - pTemp ** 2.))
    return 0.

def qManufactured1(x, p, t):
    return 0.

def infSmoothFcn(x, x1, x2):
    xTemp = (2. * x - x1 - x2) / (x2 - x1)
    return np.exp(-1. / (1. - xTemp ** 2.))

def infSmoothFcnDer(x, x1, x2):
    xTemp = (2. * x - x1 - x2) / (x2 - x1)
    return -4. * (2 * x - x1 - x2) / (x2 - x1) ** 2. \
        * 1. / (1. - xTemp ** 2.) ** 2. * np.exp(-1. / (1. - xTemp ** 2.))

def SourceManufactured1(T, q, x, p, t):
    if x >= suppLeft and x <= suppRight and p >= suppBottom and p <= suppTop:
        xTemp = (2. * x - suppLeft - suppRight) / (suppRight - suppLeft)
        pTemp = (2. * p - suppBottom - suppTop) / (suppTop - suppBottom)
        part_dTdt = -np.sin(t) \
             * np.exp(-1. / (1 - xTemp ** 2.)) \
             * np.exp(-1. / (1 - pTemp ** 2.))
        part_dudx = - np.cos(np.pi * p / p0) * 2. * np.pi / xL * np.sin(2. * np.pi * x / xL)
        part_dTdx = np.cos(t) * infSmoothFcnDer(x, suppLeft, suppRight)\
            * infSmoothFcn(p, suppBottom, suppTop) 
        part_domegadp = np.pi / p0 * np.cos(np.pi * p / p0) * np.cos(2. * np.pi * x / xL)
        part_dTdp = np.cos(t) * infSmoothFcn(x, suppLeft, suppRight)\
            * infSmoothFcnDer(p, suppBottom, suppTop) 
        return [part_dTdt + part_dudx * T + u(x, p) * part_dTdx
             + part_domegadp * T + omega(x, p) * part_dTdp, 0.]
    return [0., 0.]


# ----- Test 3 ----- #

suppLeft = x0 + 4. * (xL - x0) / 9.
suppRight = x0 + 5. * (xL - x0) / 9.
pSec = (pB(0) - pA) / 9.
suppBottom = pA + 4. * pSec
suppTop = suppBottom + pSec

def TManufactured2(x, p, t):
    if x >= suppLeft and x <= suppRight and p >= suppBottom and p <= suppTop:
        xTemp = (2. * x - suppLeft - suppRight) / (suppRight - suppLeft)
        pTemp = (2. * p - suppBottom - suppTop) / (suppTop - suppBottom)
        return np.exp(t) \
             * np.exp(-1. / (1. - xTemp ** 2.)) \
             * np.exp(-1. / (1. - pTemp ** 2.))
    return 0.

def qManufactured2(x, p, t):
    return 0.

def SourceManufactured2(T, q, x, p, t):
    if x >= suppLeft and x <= suppRight and p >= suppBottom and p <= suppTop:
        xTemp = (2. * x - suppLeft - suppRight) / (suppRight - suppLeft)
        pTemp = (2. * p - suppBottom - suppTop) / (suppTop - suppBottom)
        part_dTdt = np.exp(t) \
             * np.exp(-1. / (1 - xTemp ** 2.)) \
             * np.exp(-1. / (1 - pTemp ** 2.))
        part_dudx = - np.cos(np.pi * p / p0) * 2. * np.pi / xL * np.sin(2. * np.pi * x / xL)
        part_dTdx = np.exp(t) * infSmoothFcnDer(x, suppLeft, suppRight)\
            * infSmoothFcn(p, suppBottom, suppTop) 
        part_domegadp = np.pi / p0 * np.cos(np.pi * p / p0) * np.cos(2. * np.pi * x / xL)
        part_dTdp = np.exp(t) * infSmoothFcn(x, suppLeft, suppRight)\
            * infSmoothFcnDer(p, suppBottom, suppTop) 
        return [part_dTdt + part_dudx * T + u(x, p) * part_dTdx
             + part_domegadp * T + omega(x, p) * part_dTdp, 0.]
    return [0., 0.]