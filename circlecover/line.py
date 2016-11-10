
class line:
    """
    A class for a line segment.
    """
    def __init__(self,P1,P2):
        """
        Define a line segment.
        """
        self.P1 = P1
        self.P2 = P2
    
    def get_p1(self):
        return self.P1
        
    def get_p2(self):
        return self.P2

    def slope(self):
        """
        Returns the slope or float('inf') if line is vertical.
        """
        deltaY = float(self.P1[1] - self.P2[1])
        deltaX = float(self.P1[0] - self.P2[0])
        # A very large slope if the line is vertical.
        if deltaX == 0:
            return float('inf')
        else:
            return deltaY/deltaX

    def __str__( self ):
        return '[ ' + str(self.P1) + " , " + str(self.P2) + ' ]'
