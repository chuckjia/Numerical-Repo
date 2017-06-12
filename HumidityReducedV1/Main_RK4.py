# import numpy as np
from SpaceDiscretization import *
from Constants import *
from InitialConditions import *
import matplotlib.pyplot as plt

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Build the mesh as a 2D array of cells
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

# Initialize a mesh/a 2D array of cells

meshCells = [[Cell() for j in range(Np + 2)] for i in range(Nx + 2)]

# ----- ----- ----- ----- ----- ----- #

# Build inner cells

for i in range(1, Nx + 1):
    for j in range(1, Np + 1):

        # x coordinates
        xRight = i * Dx  # The rightmost x coordinate of cell; Can be optimized
        xLeft = xRight - Dx

        # p coordinates
        DpCurr = Dp(xRight)
        pTopRight = pA + j * DpCurr
        pBottomRight = pTopRight - DpCurr
        DpCurr = Dp(xLeft)
        pTopLeft = pA + j * DpCurr
        pBottomLeft = pTopLeft - DpCurr

        # Define the cells
        meshCells[i][j].setCell(i, j, Coord(xLeft, pTopLeft), 
                                Coord(xRight, pTopRight),
                                Coord(xLeft, pBottomLeft), 
                                Coord(xRight, pBottomRight))

# ----- ----- ----- ----- ----- ----- #

# Building the flat control volumes

# Cell (0, j)
DpCurr = Dp(x0)
for j in range(1, Np + 1):
    pTop = pA + j * DpCurr
    pBottom = pTop - DpCurr
    meshCells[0][j].setFlatCtrlVol(0, j, Coord(x0, pTop), 
                                   Coord(x0, pTop),
                                   Coord(x0, pBottom), 
                                   Coord(x0, pBottom))
# Cell (Nx + 1, j)
DpCurr = Dp(xL)
for j in range(1, Np + 1):
    pTop = pA + j * DpCurr
    pBottom = pTop - DpCurr
    meshCells[Nx + 1][j].setFlatCtrlVol(Nx + 1, j, Coord(xL, pTop), 
                                        Coord(xL, pTop),
                                        Coord(xL, pBottom), 
                                        Coord(xL, pBottom))
# Cell (i, 0)
for i in range(1, Nx + 1):
    xRight = i * Dx
    xLeft = xRight - Dx
    meshCells[i][0].setFlatCtrlVol(i, 0, Coord(xLeft, pA), 
                                   Coord(xRight, pA),
                                   Coord(xLeft, pA), 
                                   Coord(xRight, pA))
# Cell (i, Np + 1)
for i in range(1, Nx + 1):
    xRight = i * Dx
    xLeft = xRight - Dx
    pRight = pB(xRight)
    pLeft = pB(xLeft)
    meshCells[i][Np + 1].setFlatCtrlVol(i, Np + 1, Coord(xLeft, pLeft),
                                        Coord(xRight, pRight), 
                                        Coord(xLeft, pLeft),
                                        Coord(xRight, pRight))


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Initialize the solutions, flux functions, and source terms
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

# Initialize the solutions 

T = [[0. for j in range(Np + 2)] for i in range(Nx + 2)]
q = [[0. for j in range(Np + 2)] for i in range(Nx + 2)]

for i in range(1, Nx + 1):  # Initialize solutions using initial conditions
    for j in range(1, Np + 1):
        currCellCenter = meshCells[i][j].center
        xCurr = currCellCenter.x
        pCurr = currCellCenter.p
        T[i][j] = TInitial(xCurr, pCurr)
        q[i][j] = qInitial(xCurr, pCurr)

# ----- ----- ----- ----- ----- ----- #

# Initialize the fluxes G and F

GFlux1 = [[0 for j in range(Np + 1)] for i in range(Nx + 1)]
GFlux2 = [[0 for j in range(Np + 1)] for i in range(Nx + 1)]

FFlux1 = [[0 for j in range(Np + 1)] for i in range(Nx + 1)]
FFlux2 = [[0 for j in range(Np + 1)] for i in range(Nx + 1)]


# Calculation of fluxes G and F

def calcFluxes(T, q, GFlux1, GFlux2, FFlux1, FFlux2):
    """This function calculates the fluxes"""
    
    # ----- ----- ----- ----- ----- #
    # Cases: inner cells
    # ----- ----- ----- ----- ----- #
    
    for i in range(1, Nx + 1):
        for j in range(1, Np + 1):
            currCell = meshCells[i][j]
            
            # ----- ----- ----- ----- ----- #
            
            # Calculate G(i, j + 1/2): flux on the top edge of cell

            # The value uVal represents the u(i, j + 1/2) value.  
            # Method different from paper. See page 109.
            uVal = u(currCell.center.x, 
                     (currCell.topLeft.p + currCell.topRight.p) * 0.5)
            # The value omegaVal represents the omega(i, j + 1/2) value.
            omegaVal = omega(currCell.center.x, (currCell.topLeft.p 
                                                 + currCell.topRight.p) * 0.5)
            # The normal vector on the cell top as a coordinate
            normalVector = currCell.getTopNormal()
            # Product of the first 3 terms in the expression of G(i, j + 1/2)
            prodVal = currCell.getTopEdgeLength() * (normalVector.x * uVal 
                                                     + normalVector.p * omegaVal)
            # Calculate yCheck
            if prodVal >= 0:
                TCheck = T[i][j]
                qCheck = q[i][j]
            else:
                TCheck = T[i][j + 1]
                qCheck = q[i][j + 1]
            # Calculate G flux
            GFlux1[i][j] = prodVal * TCheck
            GFlux2[i][j] = prodVal * qCheck
            
            # ----- ----- ----- ----- ----- #

            # Calculating F(i + 1/2, j): flux on the right edge of cell

            xRight = i * Dx  # x_{i + 1/2}
            prodVal = Dp(xRight) * u(xRight, currCell.center.p)
            # Value of yCheck
            if prodVal >= 0:
                TCheck = T[i][j]
                qCheck = q[i][j]
            else:
                TCheck = T[i + 1][j]
                qCheck = q[i + 1][j]
            # Value of F flux
            FFlux1[i][j] = prodVal * TCheck
            FFlux2[i][j] = prodVal * qCheck

    # ----- ----- ----- ----- ----- #
    # Cases: when j = 0
    # ----- ----- ----- ----- ----- # 

    for i in range(1, Nx + 1):
        currCell = meshCells[i][0]
        # prodVal is the value of |Gamma(i + 1/2, j)|u(i + 1/2, j)
        # Note that the normal vector is (1, 0)
        prodVal = Dx * omega(currCell.center.x, pA)
        # yCheck value
        if prodVal >= 0:
            TCheck = T[i][0]
            qCheck = q[i][0]
        else:
            TCheck = T[i][1]
            qCheck = q[i][1]
        # G value
        GFlux1[i][0] = prodVal * TCheck
        GFlux2[i][0] = prodVal * qCheck
        
    # ----- ----- ----- ----- ----- #
    # Cases: when i = 0
    # ----- ----- ----- ----- ----- # 

    for j in range(1, Np + 1):
        currCell = meshCells[0][j]
        DpVal = Dp(x0)
        prodVal = DpVal * u(x0, j * DpVal)
        # Value of yCheck
        if prodVal >= 0:
            TCheck = T[0][j]
            qCheck = q[0][j]
        else:
            TCheck = T[1][j]
            qCheck = q[1][j]
        # Value of F
        FFlux1[0][j] = prodVal * TCheck
        FFlux2[0][j] = prodVal * qCheck


# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #
# Run the Runge Kutta method
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== #

# The intermediate solutions
TTemp = [[0. for j in range(Np + 2)] for i in range(Nx + 2)]
qTemp = [[0. for j in range(Np + 2)] for i in range(Nx + 2)]

# t is time variable
t = 0.
# step is the current number of step. step is for testing purposes
step = 0.

while (step < numSteps):
    # The currRunningPercentage is for testing and informational purposes
    step += 1.
    t = step * Dt
    currRunningPercentage = step / numSteps * 100
    print currRunningPercentage, "%", " Time =", t
    
    kSum1 = [[0. for j in range(Np + 2)] for i in range(Nx + 2)]
    kSum2 = [[0. for j in range(Np + 2)] for i in range(Nx + 2)]
    kCurr1 = [[0. for j in range(Np + 2)] for i in range(Nx + 2)]
    kCurr2 = [[0. for j in range(Np + 2)] for i in range(Nx + 2)] 
    
    # ----- ----- #
    # RK Step 1
    # ----- ----- #

    # Calculate G, F, S
    calcFluxes(T, q, GFlux1, GFlux2, FFlux1, FFlux2)
    
    for i in range(1, Nx + 1):
        for j in range(1, Np + 1):
            currCell = meshCells[i][j]
            x = currCell.center.x
            p = currCell.center.p
            [SVal1, SVal2] = calcS(T[i][j], q[i][j], x, p, t)
            kCurr1[i][j] = -1 / currCell.volume \
             * (GFlux1[i][j] - GFlux1[i][j - 1] + FFlux1[i][j] - FFlux1[i][j - 1]) + SVal1
            kCurr2[i][j] = -1 / currCell.volume\
             * (GFlux2[i][j] - GFlux2[i][j - 1] + FFlux2[i][j] - FFlux2[i][j - 1]) + SVal2
            kSum1[i][j] += kCurr1[i][j]
            kSum2[i][j] += kCurr2[i][j]
            TTemp[i][j] = T[i][j] + 0.5 * Dt * kCurr1[i][j]
            qTemp[i][j] = q[i][j] + 0.5 * Dt * kCurr2[i][j]

    # ----- ----- #
    # RK Step 2
    # ----- ----- #
    
    calcFluxes(TTemp, qTemp, GFlux1, GFlux2, FFlux1, FFlux2)
    
    for i in range(1, Nx + 1):
        for j in range(1, Np + 1):
            currCell = meshCells[i][j]
            x = currCell.center.x
            p = currCell.center.p
            [SVal1, SVal2] = calcS(TTemp[i][j], qTemp[i][j], x, p, t)
            kCurr1[i][j] = -1 / currCell.volume \
             * (GFlux1[i][j] - GFlux1[i][j - 1] + FFlux1[i][j] - FFlux1[i][j - 1]) + SVal1
            kCurr2[i][j] = -1 / currCell.volume\
             * (GFlux2[i][j] - GFlux2[i][j - 1] + FFlux2[i][j] - FFlux2[i][j - 1]) + SVal2
            kSum1[i][j] += 2 * kCurr1[i][j]
            kSum2[i][j] += 2 * kCurr2[i][j]
            TTemp[i][j] = T[i][j] + 0.5 * Dt * kCurr1[i][j]
            qTemp[i][j] = q[i][j] + 0.5 * Dt * kCurr2[i][j]
            
    # ----- ----- #
    # RK Step 3
    # ----- ----- #
    
    calcFluxes(TTemp, qTemp, GFlux1, GFlux2, FFlux1, FFlux2)
    
    for i in range(1, Nx + 1):
        for j in range(1, Np + 1):
            currCell = meshCells[i][j]
            x = currCell.center.x
            p = currCell.center.p
            [SVal1, SVal2] = calcS(TTemp[i][j], qTemp[i][j], x, p, t)
            kCurr1[i][j] = -1 / currCell.volume \
             * (GFlux1[i][j] - GFlux1[i][j - 1] + FFlux1[i][j] - FFlux1[i][j - 1]) + SVal1
            kCurr2[i][j] = -1 / currCell.volume\
             * (GFlux2[i][j] - GFlux2[i][j - 1] + FFlux2[i][j] - FFlux2[i][j - 1]) + SVal2
            kSum1[i][j] += 2 * kCurr1[i][j]
            kSum2[i][j] += 2 * kCurr2[i][j]
            TTemp[i][j] = T[i][j] + Dt * kCurr1[i][j]
            qTemp[i][j] = q[i][j] + Dt * kCurr2[i][j]
            
    # ----- ----- #
    # RK Step 4
    # ----- ----- #
    
    calcFluxes(TTemp, qTemp, GFlux1, GFlux2, FFlux1, FFlux2)
    
    for i in range(1, Nx + 1):
        for j in range(1, Np + 1):
            currCell = meshCells[i][j]
            x = currCell.center.x
            p = currCell.center.p
            [SVal1, SVal2] = calcS(TTemp[i][j], qTemp[i][j], x, p, t)
            kCurr1[i][j] = -1 / currCell.volume \
             * (GFlux1[i][j] - GFlux1[i][j - 1] + FFlux1[i][j] - FFlux1[i][j - 1]) + SVal1
            kCurr2[i][j] = -1 / currCell.volume\
             * (GFlux2[i][j] - GFlux2[i][j - 1] + FFlux2[i][j] - FFlux2[i][j - 1]) + SVal2
            kSum1[i][j] += kCurr1[i][j]
            kSum2[i][j] += kCurr2[i][j]
            
    for i in range(1, Nx + 1):
        for j in range(1, Np + 1):
            currCell = meshCells[i][j]
            x = currCell.center.x
            p = currCell.center.p 
            T[i][j] += 1 / 6. * Dt * kSum1[i][j]
            q[i][j] += 1 / 6. * Dt * kSum2[i][j]
    
    # ----- ----- ----- ----- ----- #
    # Boundary conditions: Neumann
    # ----- ----- ----- ----- ----- #
    
    for j in range(Np + 2):
        T[0][j] = T[1][j]
        q[0][j] = q[1][j]
        T[Nx + 1][j] = T[Nx][j]
        q[Nx + 1][j] = q[Nx][j]
        
    for i in range(Nx + 2):
        T[i][0] = T[i][1]
        q[i][0] = q[i][1]
        T[i][Np + 1] = T[i][Np]
        q[i][Np + 1] = q[i][Np]

sum1 = 0
sum2 = 0
for i in range(1, Nx + 1):
    for j in range(1, Np + 1):
        currCell = meshCells[i][j]
        x = currCell.center.x
        p = currCell.center.p
        vol = currCell.volume
        sum1 += (T[i][j] - TManufactured1(x, p, finalTime)) ** 2 * vol
        sum2 += vol * TManufactured1(x, p, finalTime) ** 2

print (sum1 / sum2) ** 0.5


plt.plot(T)       
#plt.imshow(T)
plt.show()









