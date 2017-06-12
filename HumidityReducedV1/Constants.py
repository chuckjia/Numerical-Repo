# This file stores the global constants from the humidity model, including
# functions used in the equations and in calculating the constants

import numpy as np

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# x direction constants
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

Nx = 200  # Number of cells along the x direction

x0 = 0.
xL = 50000.
Dx = xL / Nx  # Delta x


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# p direction constants
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

Np = 200  # Number of cells along the p direction

pA = 200.

def pB(x):
    """This function defines the topography at the bottom of the 
    atmosphere"""
    # See page 113 in the paper
    return 1000. - 200. * np.exp(- (x - 25000.) ** 2. / 9000000.)

def Dp(x):
    """This function calculates the Delta p"""
    return (pB(x) - pA) / Np


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Time related constants
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

numSteps = 100  # Number of time steps

Dt = 0.01  # Delta t
finalTime = numSteps * Dt  # Final time, which is used in while loops


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Other equation constants
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

p0 = 1000.  # Chosen/given in the problem setting
R = 287.  # Gas constant for dry air

# ----- ----- ----- ----- ----- #

# The function u 

def u(x, p):
    return 7.5 + np.cos(np.pi * p / p0) * np.cos(2. * np.pi * x / xL)

# ----- ----- ----- ----- ----- #

# The function omega

def omega(x, p):
    return np.sin(np.pi * p / p0) * np.cos(2. * np.pi * x / xL)


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# The source terms from the original model
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

# (This source term might not be the one used in tests)

def SourceHumidity(T, q, x, p):
    # The order of the following definitions is according to the order of
    # appearance of the variables or constants in S
    
    # Defining omega
    omegaVal = omega(x, p)
    # Defining Cp
    Cp = 1004.
    # Defining delta
    qs = 0.622 * 6.112 * np.exp(17.67 * (T - 273.15) / (T - 29.65)) / p
    delta = 0.5 * (1. - np.sign(omegaVal)) * 0.5 * (1. + np.sign(q - qs))
    # Defining L
    L = 2.5008e6 - 2.3e3 * (T - 275.)
    # Defining F
    Rv = 461.50
    F = qs * T * (L * R - Cp * Rv * T) / (Cp * Rv * T ** 2. + qs * L ** 2.)

    return [omegaVal / p * (R * T - delta * L * F) / Cp, delta / p * omegaVal * F]  

    
    
    
