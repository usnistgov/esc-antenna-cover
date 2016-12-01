from __future__ import division
import numpy as np
import math
import line
import pdb


class Circle:
    """
    A circle class.
    """
    def __init__(self, center=None, radius=None):
        """
        Constructor.
        """
        # Radius of the circle.
        self.r = radius
        # center of the rircle
        self.Q = center 

    def inside(self, point):
        """
        Return true if point is inside a circle.
        """
        try:
       	    x = point[0]
            y = point[1]
       	    return (x- self.Q[0])**2 + (y - self.Q[1])**2 <= self.r**2
        except:
            pdb.set_trace()
            

    def get_center(self):
        return self.Q

    def get_radius(self):
        return self.r

    def encloses(self, line_segment):
        """
        line_segement: A line segment that is defined by two endpoints P1 and P2
        """
        P1 = line_segment.get_p1()
        P2 = line_segment.get_p2()
        return self.inside(P1) and self.inside(P2)

    def collides(self,line_segment):
        """
        Tests if a line_segment collides with this circle.
        Return a pair (boolean,pieces).
        
        boolean : A boolean variable that indicates collision.

        if boolean is True then a collision was detected.

        if boolean is False then no collision was detected.

        - if both endpoints of the segment are inside the circle then boolean is true.
        - if both endpoints are outside, then boolean returns false.
        - if one endpoint is inside and another is outside then boolean returns True.

        pieces - a list of line segments that lie OUTSIDE the circle after
        intersection is computed. That is:

        - if the circle does not collide with the line segment, 
          pieces contains a list consisting of one element which is the 
          original line segment.
        - if the circle includes both endpoints then pieces is an empty list.
        - if the circle includes only one endpoint then pieces is a list of one line segment.

        """
        

        if self.encloses(line_segment):
            # Circle encloses line.  Return empty list.
            return True,[]

        # Translate the center of the circle to 0
        # The line segment moves as a result.
        x1 = line_segment.get_p1()[0] - self.Q[0]
        x2 = line_segment.get_p2()[0] - self.Q[0]
        y1 = line_segment.get_p1()[1] - self.Q[1]
        y2 = line_segment.get_p2()[1] - self.Q[1]
        dx = x2 - x1
        dy = y2 - y1
        dr_sqr = dx**2 + dy**2
        D = x1*y2 - x2*y1
        # the discrimininant 
        disc = (self.r**2)*(dr_sqr) - D**2
        # discriminant is < 0 then no intersection.
        # discriminant is 0 then tangent.
        if disc < 0 or disc == 0:
            return False,[line_segment]
        sqrt_disc = math.sqrt(disc)
        sign_dy = 1
        if dy < 0:
            sign_dy = -1
        # translate the points back to original locations. 
        ix0 = (D*dy + sign_dy*dx*sqrt_disc)/dr_sqr + self.Q[0]
        ix1 = (D*dy - sign_dy*dx*sqrt_disc)/dr_sqr + self.Q[0]
        iy0 = (-D*dx + abs(dy)*sqrt_disc)/dr_sqr + self.Q[1]
        iy1 = (-D*dx - abs(dy)*sqrt_disc)/dr_sqr + self.Q[1]
        # the intersection points are known but we need 
        # to order them.
        pointlist = []
        pointlist.append(line_segment.get_p1())
        pointlist.append([ix0,iy0])
        pointlist.append([ix1,iy1])
        pointlist.append(line_segment.get_p2())
        # sort in the same order as the line segment:
        if line_segment.get_sort_dimension() == 0:
            pointlist.sort(key=lambda x : x[0])
        else:
            pointlist.sort(key=lambda x : x[1])
        # We now have 4 colinear points. 
        # Check the location of the original points of the line segment.
        retval = []
        res = True
        if line_segment.get_p1() == pointlist[0] and \
            line_segment.get_p2() == pointlist[1] :
            # The line segment would intersect with circle if extended
            res = False
	    retval.append(line_segment)
        elif line_segment.get_p1() == pointlist[0] and \
            line_segment.get_p2() == pointlist[3] :
            # In this case there are 2 segments remaining
	    if not np.allclose(pointlist[0],pointlist[1]):
            	l1 = line.Line(pointlist[0],pointlist[1])
            	retval.append(l1)
	    if not np.allclose(pointlist[2],pointlist[3]):
            	l2 = line.Line(pointlist[2],pointlist[3])
            	retval.append(l2)
        elif line_segment.get_p1() == pointlist[0] :
	    if not np.allclose(pointlist[0],pointlist[1]):
            	l1 = line.Line(pointlist[0],pointlist[1])
            	retval.append(l1)
        elif line_segment.get_p2() == pointlist[3] :
	    if not np.allclose(pointlist[2],pointlist[3]):
            	l2 = line.Line(pointlist[2],pointlist[3])
            	retval.append(l2)
        return res,retval



    def __hash__(self):
        return hash(str(self.Q) + str(self.r))

    def __eq__(self,other):
        return self.Q == other.Q and self.r == other.r

    def __repr__(self):
        return '{ center : ' + str(self.Q) + ', radius: ' + str(self.r) + '}'
