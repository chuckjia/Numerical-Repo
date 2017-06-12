class Coord(object):
    """This class represents the coordinates of a grid point"""

    def __init__(self, x, p):
        self.x = x
        self.p = p

    def distFrom(self, coord):
        return ((self.x - coord.x) ** 2 + (self.p - coord.p) ** 2) ** 0.5

    def __str__(self):
        return "Coordinate (" + str(self.x) + ", " + str(self.p) + "): x = " + str(self.x) + ", p = " + str(self.p)


class Cell(object):
    """This class represents a cell. i and j are indices of the cell.
    topLeft, topRight, bottomLeft, bottomRight are Coord"""

    def __init__(self):
        self.i = None
        self.j = None
        self.topLeft = None
        self.topRight = None
        self.bottomLeft = None
        self.bottomRight = None
        self.center = None
        self.volume = None
        
    # Setter functions

    def setIndex(self, i, j):
        self.i = i
        self.j = j

    def setGeometry(self, topLeft, topRight, bottomLeft, bottomRight):
        self.topLeft = topLeft
        self.topRight = topRight
        self.bottomLeft = bottomLeft
        self.bottomRight = bottomRight
        self.volume = self.setVolume()

    def setBaryCenter(self):
        """This method calculates the barycenter of a cell using the method given in
        the paper"""
        # South-East triangle
        xbar1 = (self.bottomLeft.x + self.bottomRight.x + self.topRight.x) / 3.0
        pbar1 = (self.bottomLeft.p + self.bottomRight.p + self.topRight.p) / 3.0
        # North-East triangle
        xbar3 = (self.topLeft.x + self.bottomRight.x + self.topRight.x) / 3.0
        pbar3 = (self.topLeft.p + self.bottomRight.p + self.topRight.p) / 3.0
        # North-West triangle
        xbar2 = (self.topLeft.x + self.bottomLeft.x + self.topRight.x) / 3.0
        pbar2 = (self.topLeft.p + self.bottomLeft.p + self.topRight.p) / 3.0
        # South-West triangle
        xbar4 = (self.topLeft.x + self.bottomLeft.x + self.bottomRight.x) / 3.0
        pbar4 = (self.topLeft.p + self.bottomLeft.p + self.bottomRight.p) / 3.0
        # Calculate the intersection of two diagonal lines
        denom = (xbar1 - xbar2) * (pbar3 - pbar4) - (pbar1 - pbar2) * (xbar3 - xbar4)
        x = ((xbar1 * pbar2 - pbar1 * xbar2) * (xbar3 - xbar4) - (xbar1 - xbar2) * (
            xbar3 * pbar4 - pbar3 * xbar4)) / denom
        y = ((xbar1 * pbar2 - pbar1 * xbar2) * (pbar3 - pbar4) - (pbar1 - pbar2) * (
            xbar3 * pbar4 - pbar3 * xbar4)) / denom
        self.center = Coord(x, y)

    def setCenter(self):
        """This method calculates the center of a flat control volume"""
        x = (self.topLeft.x + self.topRight.x + self.bottomLeft.x + self.bottomRight.x) / 4.0
        y = (self.topLeft.p + self.topRight.p + self.bottomLeft.p + self.bottomRight.p) / 4.0
        self.center = Coord(x, y)

    def setVolume(self):
        """This method calculates the size (volume or area) of the cell"""
        return (self.topRight.x - self.topLeft.x) * (
            self.topLeft.p - self.bottomLeft.p + self.topRight.p - self.bottomRight.p) / 2.0

    def setCell(self, i, j, topLeft, topRight, bottomLeft, bottomRight):
        """Set all the cell properties for an inner (normal) cell (meaning 
           not a flat control volume)"""
        self.setIndex(i, j)
        self.setGeometry(topLeft, topRight, bottomLeft, bottomRight)
        self.setBaryCenter()

    def setFlatCtrlVol(self, i, j, topLeft, topRight, bottomLeft, bottomRight):
        """Set all the cell properties for a flat control volume)"""
        self.setIndex(i, j)
        self.setGeometry(topLeft, topRight, bottomLeft, bottomRight)
        self.setCenter()

    # Getter functions

    def getTopNormal(self): # Might have problems on the direction
        x = self.topRight.x - self.topLeft.x
        p = self.topRight.p - self.topLeft.p
        denom = (x ** 2 + p ** 2) ** 0.5
        return Coord(-p / denom, x / denom)

    def getLeftNormal(self):
        return Coord(1, 0)

    def getTopEdgeLength(self):
        return ((self.topRight.x - self.topLeft.x) ** 2
                + (self.topRight.p - self.topLeft.p) ** 2) ** 0.5

    def getBottomEdgeLength(self):
        return ((self.bottomRight.x - self.bottomLeft.x) ** 2
                + (self.bottomRight.p - self.bottomLeft.p) ** 2) ** 0.5

    def __str__(self):
        return "Cell (%d, %d)" % (self.i, self.j)
