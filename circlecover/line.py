
# This software was developed by employees of the National Institute
# of Standards and Technology (NIST), an agency of the Federal
# Government. Pursuant to title 17 United States Code Section 105, works
# of NIST employees are not subject to copyright protection in the United
# States and are considered to be in the public domain. Permission to freely
# use, copy, modify, and distribute this software and its documentation
# without fee is hereby granted, provided that this notice and disclaimer
# of warranty appears in all copies.
#
# THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
# EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
# TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
# IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
# AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION
# WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE
# ERROR FREE. IN NO EVENT SHALL NASA BE LIABLE FOR ANY DAMAGES, INCLUDING,
# BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES,
# ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS
# SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR
# OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY
# OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT
# OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
#
# Distributions of NIST software should also include copyright and licensing
# statements of any third-party software that are legally bundled with
# the code in compliance with the conditions of those licenses.

import math
import numpy as np
import pdb
from shapely.geometry import LineString
from shapely.geometry import MultiPoint
import json

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


    def intersects(self,line2):
        """
        Fast boolean check to test if given line segment
        intersects another line segment.
        http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
        """
        def ccw(A,B,C):
            return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

        # Return true if line segments AB and CD intersect
        A = self.coords[0] 
        B = self.coords[1]
        C = line2.coords[0]
        D = line2.coords[1]
        return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)    


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

    def get_centroid(self):
        return ((self.coords[0][0] + self.coords[1][0])/2 , (self.coords[0][1] + self.coords[1][1])/2)

    def find_centroid_of_closest_lines(self,lines,count):
        my_centroid = self.get_centroid()
        # Sort the distances from my centroid to the line set
        sorted_distances = sorted([(k,l,l.distance_to_centroid(my_centroid)) for (k,l) in enumerate(lines) if l != self ], key = lambda t: t[2])
        if len(sorted_distances) > count :
            s =  sorted_distances[0:count]
        else:
            s = sorted_distances
        # Return my k nearest neighbors
        points = [(l[1].get_p1(),l[1].get_p2()) for  l in s ]
        mp = MultiPoint(points)
        return mp.centroid.coords[0]
        


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
            

    def get_coordinates(self):
        return [ self.coords[0],self.coords[1] ]


    def length(self):
        x1 = self.coords[0][0]
        y1 = self.coords[0][1]
        x2 = self.coords[1][0]
        y2 = self.coords[1][1]
        return math.sqrt((x1 - x2)**2 + (y1-y2)**2)


    def distance(self,point):
        """
        distance between a point and a line.
        """
        y1 = self.get_p1()[1]
        y2 = self.get_p2()[1]
        x1 = self.get_p1()[0]
        x2 = self.get_p2()[1]
        x0 = point[0]
        y0 = point[1]
        return abs((y2 -y1)*x0 - (x2 -x1)*y1 + x2*y1 -y2*x1)/ math.sqrt((y2-y1)**2 + (x2-x1)**2)

    
    def circumscribing_radius(self,point):
        return math.sqrt(max( 
                ((point[0]  - self.get_p1()[0])**2 + (point[1] - self.get_p1()[1])**2),
                ((point[0]  - self.get_p2()[0])**2 + (point[1] - self.get_p2()[1])**2)))
    
    def distance_to_centroid(self,point):
        mp = [(self.get_p1()[0] + self.get_p2()[0])/2.0,(self.get_p1()[1] + self.get_p2()[1])/2.0]
        return math.sqrt((point[0] - mp[0])**2 + (point[1] - mp[1])**2)


    def __hash__(self):
        return hash(str(self.coords))


    def __repr__( self ):
        return "Line : " + json.dumps([[self.coords[0][0],self.coords[0][1]],[self.coords[1][0],self.coords[1][1]]]) 

    def __eq__(self,other):
        if other == None: 
            return False
        return np.allclose(self.get_p1(),other.get_p1(),atol=.001) and np.allclose(self.get_p2(),other.get_p2(),atol=.001)
