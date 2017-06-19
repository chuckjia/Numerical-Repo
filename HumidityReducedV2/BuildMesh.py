from MeshFunctions import *

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== # 
# Build the Mesh
# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== # 

# ----- ----- ----- ----- #
# Build the Grid
# ----- ----- ----- ----- #

# Initialize 2D array of y coord for the top right corner of cell[i][j]
meshGridY = np.zeros((Nx + 1, Np + 1)) 

# Calculate grid y coord
# For the boundary i = 0    
DpVal = Dp(x0)
for j in range(1, Np + 1):
    meshGridY[0][j] = calcCellTopRightY(j, DpVal)
# For the boundary j = 0  
DpVal = Dp(xL)
for i in range(0, Nx + 1):
    meshGridY[i][0] = pA
        
# Build the mesh for all other cells
for i in range(1, Nx + 1):
    xVal = i * Dx
    DpVal = Dp(xVal)
    for j in range(1, Np + 1):
        meshGridY[i][j] = calcCellTopRightY(j, DpVal)
    
        
def getCellTopRight(i, j):
    if j == Np + 1:
        j -= 1
    return [calcCellRightX(i), meshGridY[i][j]]

def getCellTopLeft(i, j):
    if i == 0:
        i += 1
    return [calcCellLeftX(i), meshGridY[i - 1][j]]

def getCellBottomRight(i, j):
    if j == 0:
        j += 1
    return getCellTopRight(i, j - 1);

def getCellBottomLeft(i, j):
    if i == 0:
        i += 1
    if j == 0:
        j += 1
    return getCellTopRight(i - 1, j - 1);
  
# ----- ----- ----- ----- ----- ----- #
# Calculate Volumes and Bary Centers
# ----- ----- ----- ----- ----- ----- #

# Initialize 2D array of volumes for the cell[i][j]
meshCellVol = np.zeros((Nx + 2, Np + 2))
meshCellCenterX = np.zeros((Nx + 2, Np + 2))
meshCellCenterY = np.zeros((Nx + 2, Np + 2))

for i in range(Nx + 2):
    for j in range(Np + 2):
        topRight = getCellTopRight(i, j)
        topLeft = getCellTopLeft(i, j)
        bottomLeft = getCellBottomLeft(i, j)
        bottomRight = getCellBottomLeft(i, j)
        meshCellVol = calcVol(topLeft, topRight, bottomLeft, bottomRight)
        cellCenter = calcBaryCenter(topLeft, topRight, bottomLeft, bottomRight)
        meshCellCenterX[i][j] = cellCenter[0]
        meshCellCenterY[i][j] = cellCenter[1]
        