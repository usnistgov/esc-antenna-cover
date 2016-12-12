import math
import numpy as np
import pdb
from shapely.geometry import LineString

class Line(LineString):
    """
    A class for a line segment.
    """
    def __init__(self,p1,p2):
        """
        Define a line segment.
        """
        points = [p1,p2]
        if p1[0] != p2[0]:
            points.sort(key=lambda x : x[0])
            self.sort_dimension = 0
        else:
            points.sort(key=lambda x : x[1])
            self.sort_dimension = 1
        self.circles = []
        LineString.__init__(self,points)

    def get_sort_dimension(self):
        return self.sort_dimension
    
    def get_p1(self):
        return self.coords[0]
        
    def get_p2(self):
        return self.coords[1]

    def slope(self):
        """
        Returns the slope or float('inf') if line is vertical.
        """
        deltaY = float(self.coords[0][1] - self.coords[1][1])
        deltaX = float(self.coords[0][0] - self.coords[1][0])
        # A very large slope if the line is vertical.
        if deltaX == 0:
            return float('inf')
        else:
            return deltaY/deltaX

    def isCollinear(self,point):
        """
        return True if the given point is 
        collinear with this line.
        """
        l1 = Line(self.coords[0],point)
        l2 = Line(self.coords[1],point)
        a = l1.length()
        b = l2.length()
        c = self.length()
        s = (a + b + c) /2
        return np.allclose( s-a,0) or np.allclose(s-b,0) or np.allclose(s-c, 0)

    def split(self, segmentCount):
        """
        Split a line into multiple segments.
        segmentCount : # of segments to split the line into.
        """
        segment_list = []
        deltaY = float(self.coords[1][1] - self.coords[0][1])/float(segmentCount)
        deltaX = float(self.coords[1][0] - self.coords[0][0])/float(segmentCount)
        pi_0 = self.coords[0]
        for i in range(1,segmentCount):
            yi = self.coords[0][1] + deltaY*i
            xi = self.coords[0][0] + deltaX*i
            pi_1 = [xi,yi]
            segment_list.append(Line(pi_0,pi_1))
            pi_0 = pi_1
        pi_1 = self.coords[1]
        segment_list.append(Line(pi_0,pi_1))
        return segment_list

    
    def angle(self,lineB):
        """
        return the angle between the given line segment and another line segment.

        http://stackoverflow.com/questions/28260962/calculating-angles-between-line-segments-python-with-math-atan2

        """
        def dot(vA, vB):
            return vA[0]*vB[0]+vA[1]*vB[1]
        # Dx, Dy representation for the two lines.
        vA = [(self.coords[0][0]-self.coords[1][0]), (self.coords[0][1]-self.coords[1][1])]
        vB = [(lineB.coords[0][0]-lineB.coords[1][0]), (lineB.coords[0][1]-lineB.coords[1][1])]
        # Get dot prod
        dot_prod = dot(vA, vB)
        # Get magnitudes
        magA = dot(vA, vA)**0.5
        magB = dot(vB, vB)**0.5
        # Get angle in radians 
        angle = math.acos(dot_prod/(magB*magA))
        if angle > math.pi:
            return 2*math.pi - angle
        else: 
            return angle

    def add_enclosing(self,circ):
        self.circles.append(circ)

    def clear_enclosing(self):
        self.circles = []
        

    def intersection(self, line2):
        """
        Return True and the intersection point if this line
        intersects with the given line segment line2 
        """
        if not self.intersects(line2):
            return False,None
        x = LineString.intersection(self,line2)
        return True,x.coords[0]
            

    def get_coordinates():
        return [ self.coords[0],self.coords[1] ]

    def length(self):
        x1 = self.coords[0][0]
        y1 = self.coords[0][1]
        x2 = self.coords[1][0]
        y2 = self.coords[1][1]
        return math.sqrt((x1 - x2)**2 + (y1-y2)**2)

    def __hash__(self):
        return hash(str(self.coords))

    def __repr__( self ):
        return "Line: " + str([[self.coords[0][0],self.coords[0][1]],[self.coords[1][0],self.coords[1][1]]]) + "\n"

    def __eq__(self,other):
        if other == None: 
            return False
        return np.allclose(self.coords[0],other.coords[0],atol=.001) and np.allclose(self.coords[1],other.coords[1],atol=.001)
