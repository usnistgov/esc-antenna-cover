import math
import numpy as np
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

    
    def angle(self,lineB):
        """
        return the angle between the given line and another line.

        http://stackoverflow.com/questions/28260962/calculating-angles-between-line-segments-python-with-math-atan2

        """
        def dot(vA, vB):
            return vA[0]*vB[0]+vA[1]*vB[1]
        # Get nicer vector form
        vA = [(self.points[0][0]-self.points[1][0]), (self.points[0][1]-self.points[1][1])]
        vB = [(lineB.points[0][0]-lineB.points[1][0]), (lineB.points[0][1]-lineB.points[1][1])]
        # Get dot prod
        dot_prod = dot(vA, vB)
        # Get magnitudes
        magA = dot(vA, vA)**0.5
        magB = dot(vB, vB)**0.5
        # Get angle in radians 
        angle = math.acos(dot_prod/(magB*magA))
        return angle

        
            

    def intersection(self, line2):
        """
        Return True and the intersection point if this line
        intersects with the given line segment line2
        
        Derived from:

        http://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python
        """
        xdiff = (self.points[0][0] - self.points[1][0], line2.points[0][0] - line2.points[1][0])
        ydiff = (self.points[0][1] - self.points[1][1], line2.points[0][1] - line2.points[1][1])

        def det(a, b):
            return a[0] * b[1] - a[1] * b[0]

        div = det(xdiff, ydiff)
        if div == 0:
            return False, None

        d = (det(*tuple(self.points)), det(*tuple(line2.points)))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div

        # Intersection point should be between the endpoints of the 
        # two segments.
        pointlist = []
        pointlist.append(line2.get_p1())
        pointlist.append([x,y])
        pointlist.append(line2.get_p2())
        # sort in the same order as the line segment:
        if line2.get_sort_dimension() == 0:
            pointlist.sort(key=lambda x : x[0])
        else:
            pointlist.sort(key=lambda x : x[1])

        if pointlist[0] != line2.get_p1() or \
            pointlist[2] != line2.get_p2():
            return False,None
            
        pointlist = []
        pointlist.append(self.get_p1())
        pointlist.append([x,y])
        pointlist.append(self.get_p2())
        # sort in the same order as the line segment:
        if self.get_sort_dimension() == 0:
            pointlist.sort(key=lambda x : x[0])
        else:
            pointlist.sort(key=lambda x : x[1])

        if pointlist[0] != self.get_p1() or \
            pointlist[2] != self.get_p2():
            return False,None

        return True,[x,y]


    def length(self):
        x1 = self.points[0][0]
        y1 = self.points[0][1]
        x2 = self.points[1][0]
        y2 = self.points[1][1]
        return math.sqrt((x1 - x2)**2 + (y1-y2)**2)

    def __hash__(self):
        return hash(str(self.points))

    def __repr__( self ):
        return str(self.points)

    def __eq__(self,other):
        if other == None: 
            return False
        return np.allclose(self.points[0],other.points[0],atol=.01) and np.allclose(self.points[1],other.points[1],atol=.01)
