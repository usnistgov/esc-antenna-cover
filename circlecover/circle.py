from __future__ import division
import numpy as np
import numpy.linalg.linalg as la
import math
import line


class Circle:
    """
    A circle class.
    See http://codereview.stackexchange.com/questions/86421/line-segment-to-circle-collision-algorithm    
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
        x = point[0]
        y = point[1]
        return (x- self.Q[0])**2 + (y - self.Q[1])**2 <= self.r**2

    def encloses(self, line_segment):
        """
        line_segemnt: A line segment that is defined by two endpoints P1 and P2
        """
        P1 = line_segment.get_p1()
        P2 = line_segment.get_p2()
        return self.inside(P1) and self.inside(P2)

    def collides(self,line_segment):

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
            return False,None
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
        # sorty by x and then by y
        pointlist.sort(key=lambda x : x[0])
        pointlist.sort(key=lambda x : x[1])
        # We now have 4 colinear points. 
        # Check the location of the original points of the line segment.
        retval = []
        res = True
        if line_segment.get_p1() == pointlist[0] and \
            line_segment.get_p2() == pointlist[1] :
            # The line segment would intersect with circle if extended
            res = False
        elif line_segment.get_p1() == pointlist[0] and \
            line_segment.get_p2() == pointlist[3] :
            # In this case there are 2 segments remaining
            l1 = line.Line(pointlist[0],pointlist[1])
            l2 = line.Line(pointlist[2],pointlist[3])
            retval = []
            retval.append(l1)
            retval.append(l2)
        elif line_segment.get_p1() == pointlist[0] :
            l1 = line.Line(pointlist[0],pointlist[1])
            retval = []
            retval.append(l1)
        elif line_segment.get_p2() == pointlist[3] :
            l2 = line.Line(pointlist[2],pointlist[3])
            retval.append(l2)
        return res,retval


    def collides1(self, line_segment):
        """
        line_segemnt: A line segment that is defined by two endpoints P1 and P2

        
        This returns the list of segments OUTSIDE the circle and a boolean which indicates
        whether or not the circle collides with or encloses the line segment.

        """
        if self.encloses(line_segment):
            # Circle encloses line.  Return empty list.
            return True,[]
        P1 = line_segment.get_p1()
        P2 = line_segment.get_p2()
        V = np.subtract(P2,P1)
        a = la.dot(V,V)
        Q = self.Q
        r = self.r
        b = 2 * la.dot(V,np.subtract(P1,self.Q))
        c = la.dot(P1,P1) + la.dot(Q,Q) - 2 * la.dot(P1,Q) - r**2
        disc = b**2 - 4 * a * c
        #If discriminant is negative, then there are no real solutions to the quadratic equation;
        #that means that the line misses the circle entirely.
        if disc < 0:
            return False, None
        sqrt_disc = math.sqrt(disc)
        t1 = (-b + sqrt_disc) / (2 * a)
        t2 = (-b - sqrt_disc) / (2 * a)
        if not (0 <= t1 <= 1)  and not(0 <= t2 <= 1):
            #If neither discriminant is between 0 and 1, then the line segment misses the 
            #circle (but would hit it if extended)        
            return False, None
        elif t1 <= 1 and t1 >= 0 and not (0 <= t2 <= 1):
            # Check the endpoint inside the line.
            if self.inside(P1):
                p = P2
            else:
                p = P1
            return True, [line.Line(p,np.add(P1,np.multiply(t1 , V)).tolist())]
        elif t2<=1 and t2>=0 and not (0 <= t1 <= 1):
            if self.inside(P1):
                p = P2
            else:
                p = P1
            return True, [line.Line(p,np.add(P1,np.multiply(t2 , V)).tolist())]
        else:
            q1 =  line.Line(P1,np.add(P1,np.multiply(t2 , V)).tolist())
            q2 =  line.Line(P2,np.add(P1,np.multiply(t1 , V)).tolist())
            return True, [q1,q2]

    def __repr__(self):
        return '{ center : ' + str(self.Q) + ', radius: ' + str(self.r) + '}'
