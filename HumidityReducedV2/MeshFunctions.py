from Constants import *
from Tkconstants import BOTTOM

def calcCellRightX(i):
    if i == Nx + 1:
        return xL;
    return x0 + Dx * i

def calcCellLeftX(i):
    if i == 0:
        return x0;
    return x0 + Dx * (i - 1)

def calcCellTopRightY(i, DpVal):
    if i == Np + 1:
        return pA + Np * DpVal
    return pA + DpVal * i

def calcBaryCenter(topLeft, topRight, bottomLeft, bottomRight):
    # South-East triangle
    xbar1 = (bottomLeft[0] + bottomRight[0] + topRight[0]) / 3.0
    pbar1 = (bottomLeft[1] + bottomRight[1] + topRight[1]) / 3.0
    # North-East triangle
    xbar3 = (topLeft[0] + bottomRight[0] + topRight[0]) / 3.0
    pbar3 = (topLeft[1] + bottomRight[1] + topRight[1]) / 3.0
    # North-West triangle
    xbar2 = (topLeft[0] + bottomLeft[0] + topRight[0]) / 3.0
    pbar2 = (topLeft[1] + bottomLeft[1] + topRight[1]) / 3.0
    # South-West triangle
    xbar4 = (topLeft[0] + bottomLeft[0] + bottomRight[0]) / 3.0
    pbar4 = (topLeft[1] + bottomLeft[1] + bottomRight[1]) / 3.0
    # Calculate the intersection of two diagonal lines
    denom = (xbar1 - xbar2) * (pbar3 - pbar4) - (pbar1 - pbar2) * (xbar3 - xbar4)
    x = ((xbar1 * pbar2 - pbar1 * xbar2) * (xbar3 - xbar4) - (xbar1 - xbar2) * (
        xbar3 * pbar4 - pbar3 * xbar4)) / denom
    y = ((xbar1 * pbar2 - pbar1 * xbar2) * (pbar3 - pbar4) - (pbar1 - pbar2) * (
        xbar3 * pbar4 - pbar3 * xbar4)) / denom
    return [x, y]

def calcVol(topLeft, topRight, bottomLeft, bottomRight):
    return 0.5 * (topLeft - bottomLeft + topRight - bottomRight) * (topRight[0] - topLeft[0])
