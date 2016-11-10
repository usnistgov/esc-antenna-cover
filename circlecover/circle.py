from __future__ import division
import numpy as np
import numpy.linalg.linalg as la
import math
import line


class circle:
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
        

    def collides(self, line_segment):
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
            return True, [line.line(p,np.add(P1,np.multiply(t1 , V)).tolist())]
        elif t2<=1 and t2>=0 and not (0 <= t1 <= 1):
            if self.inside(P1):
                p = P2
            else:
                p = P1
            return True, [line.line(p,np.add(P1,np.multiply(t2 , V)).tolist())]
        else:
            q1 =  line.line(P1,np.add(P1,np.multiply(t2 , V)).tolist())
            q2 =  line.line(P2,np.add(P1,np.multiply(t1 , V)).tolist())
            return True, [q1,q2]
