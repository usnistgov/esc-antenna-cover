import math
class Line:
    """
    A class for a line segment.
    """
    def __init__(self,P1,P2):
        """
        Define a line segment.
        """
        self.points =  []
        self.points.append(P1)
        self.points.append(P2)
        # sort the line segment coordinates
        if self.points[0][0] != self.points[1][0]:
       	    self.points.sort(key=lambda x : x[0])
            self.sort_dimension = 0
        else:
            self.points.sort(key=lambda x : x[1])
            self.sort_dimension = 1
    
    def get_sort_dimension(self):
        return self.sort_dimension

    def get_p1(self):
        return self.points[0]
        
    def get_p2(self):
        return self.points[1]

    def slope(self):
        """
        Returns the slope or float('inf') if line is vertical.
        """
        deltaY = float(self.points[0][1] - self.points[1][1])
        deltaX = float(self.points[0][0] - self.points[1][0])
        # A very large slope if the line is vertical.
        if deltaX == 0:
            return float('inf')
        else:
            return deltaY/deltaX

    def split(self, segmentCount):
        """
        Split a line into multiple segments.
	segmentCount : # of segments to split the line into.
        """
        segment_list = []
        deltaY = float(self.points[1][1] - self.points[0][1])/float(segmentCount)
        deltaX = float(self.points[1][0] - self.points[0][0])/float(segmentCount)
        pi_0 = self.points[0]
        for i in range(1,segmentCount):
            yi = self.points[0][1] + deltaY*i
            xi = self.points[0][0] + deltaX*i
            pi_1 = [xi,yi]
            segment_list.append(Line(pi_0,pi_1))
            pi_0 = pi_1
        pi_1 = self.points[1]
        segment_list.append(Line(pi_0,pi_1))
        return segment_list
            

    def length(self):
        x1 = self.points[0][0]
        y1 = self.points[1][1]
        x2 = self.points[1][0]
        y2 = self.points[1][1]
        return sqrt((x1 - x2)**2 + (y1-y2)**2)

    def __repr__( self ):
        return str(self.points)
