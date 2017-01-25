from __future__ import division
import numpy as np
import math
import line
import pdb
import random
import json
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import MultiLineString




class Circle:
    """
    A circle class. This is an extension of the
    shapely point class but we implement a few methods
    that are specific to circles.
    """
    def __init__(self, center=None, radius=None):
        """
        Constructor.
        """
        self.point = Point(center)
        # Radius of the circle.
        self.r = radius
        # center of the rircle
        # store this as a tuple.
        self.center = tuple(center)

    def intersects_line(self, other):
        return self.get_geometry().intersects(other)

    def distance_from_boundary(self,line):
        x1 = line.get_p1()[0]
        y1 = line.get_p1()[1]
        x2 = line.get_p2()[0]
        y2 = line.get_p2()[1]
        x,y = (x1+x2)/2 , (y1+y2)/2
        d = math.sqrt((x- self.center[0])**2 + (y - self.center[1])**2 )
        # Note that this can return a negative
        # negative means outside
        return d - self.r

    def overlaps(self,other):
        x1 = self.get_center()[0]
        y1 = self.get_center()[1]
        x2 = other.get_center()[0]
        y2 = other.get_center()[1]
        return (x2 -x1) ** 2 + (y2-y1)**2 < (self.get_radius() + other.get_radius())**2
        

    def inside(self, point):
        """
        Return true if point is inside a circle.
        """
        x = point[0]
        y = point[1]
        return (x- self.center[0])**2 + (y - self.center[1])**2 <=  1.0005*self.r**2  

    def set_radius(self, newradius):
         self.r = newradius

    def get_geometry(self):
        return self.point.buffer(self.r)

    def on(self,point):
        """
        Return True if the point is on the circle.
        """
        dist = math.sqrt((point[0] - self.get_center()[0])**2 + (point[1] - self.get_center()[1])**2)
        return np.allclose(dist,self.r,atol=.0001)
            

    def area(self):
        return math.pi*self.r**2

    def get_center(self):
        return self.center

    def get_radius(self):
        return self.r

    def encloses(self, line_segment):
        """
        line_segement: A line segment that is defined by two endpoints P1 and P2
        """
        p1 = line_segment.get_p1()
        p2 = line_segment.get_p2()
        return self.inside(p1) and self.inside(p2)

    def collides(self,line_segment):
        """
        Tests if a line_segment collides with this circle.
        Return a pair (boolean,pieces,covered).
        
        boolean : A boolean variable that indicates collision.

        if boolean is True then a collision was detected.

        if boolean is False then no collision was detected.

        - if both endpoints of the segment are inside the circle then boolean is true.
        - if both endpoints are outside, then boolean returns false.
        - if one endpoint is inside and another is outside then boolean returns True.

        pieces - a list of line segments that lie OUTSIDE the circle after
        intersection is computed. pieces union covered is line_segment
        segment. That is:
    
        covered - a line segment that is covered by the circle.

        - if the circle does not collide with the line segment, 
          "pieces" contains a list consisting of one element which is the 
          original line segment and "covered" is None.
        - if the circle includes both endpoints then pieces is an empty list.
        - if the circle includes only one endpoint then pieces is a list of one line segment 
          and "covered" is the rest of the line..

        """
        

        if self.encloses(line_segment):
            # Circle encloses line.  Return empty list.
            return True,[],line_segment

        # Translate the center of the circle to 0
        # The line segment moves as a result.
        x1 = line_segment.get_p1()[0] - self.center[0]
        x2 = line_segment.get_p2()[0] - self.center[0]
        y1 = line_segment.get_p1()[1] - self.center[1]
        y2 = line_segment.get_p2()[1] - self.center[1]
        dx = x2 - x1
        dy = y2 - y1
        dr_sqr = dx**2 + dy**2
        D = x1*y2 - x2*y1
        # the discrimininant 
        disc = (self.r**2)*(dr_sqr) - D**2
        # discriminant is < 0 then no intersection.
        # discriminant is 0 then tangent.
        # In both cases, there is no part of the line that
        # is covered.
        if disc < 0 or disc == 0:
            return False,[line_segment],None
        sqrt_disc = math.sqrt(disc)
        sign_dy = 1
        if dy < 0:
            sign_dy = -1
        # ix and iy are the coordinates of the intersection
        # of the circle with the extended line segment.
        # translate the points back to original locations. 
        ix0 = (D*dy + sign_dy*dx*sqrt_disc)/dr_sqr + self.center[0]
        ix1 = (D*dy - sign_dy*dx*sqrt_disc)/dr_sqr + self.center[0]
        iy0 = (-D*dx + abs(dy)*sqrt_disc)/dr_sqr + self.center[1]
        iy1 = (-D*dx - abs(dy)*sqrt_disc)/dr_sqr + self.center[1]
        # the intersection points are now known but we need 
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
        excluded = []
        # covered is a set of segments inside the circle.
        covered = None
        res = True
        if line_segment.get_p1() == pointlist[0] and \
            line_segment.get_p2() == pointlist[1] :
            # The line segment would intersect with circle if extended
            # but it is not covered by the circle. So the covered set is null.
            res = False
            excluded.append(line_segment)
        elif line_segment.get_p1() == pointlist[2] and \
            line_segment.get_p2() == pointlist[3] :
            # The line segment would intersect with circle if extended
            # but it is not covered by the circle. So the covered set is null.
            res = False
            excluded.append(line_segment)
        elif line_segment.get_p1() == pointlist[0] and \
            line_segment.get_p2() == pointlist[3] :
            # In this case there are 2 segments remaining
            # one segment has been covered by the line.
            if not np.allclose(pointlist[0],pointlist[1],rtol=.005):
                l1 = line.Line(pointlist[0],pointlist[1])
                excluded.append(l1)
            if not np.allclose(pointlist[2],pointlist[3],rtol=.005):
                l2 = line.Line(pointlist[2],pointlist[3])
                excluded.append(l2)
            if True or not np.allclose(pointlist[1],pointlist[2],rtol=.005):
                covered = line.Line(pointlist[1],pointlist[2])
        elif line_segment.get_p1() == pointlist[0] and \
             line_segment.get_p2() == pointlist[2] :
            # there is one segment outside the circle.
            # the rest of the line is inside the ecircle.
            l1 = line.Line(pointlist[0],pointlist[1])
            excluded.append(l1)
            if True or not np.allclose(pointlist[1],pointlist[2]):
                covered = line.Line(pointlist[1],pointlist[2])
        elif line_segment.get_p1() == pointlist[1] and \
            line_segment.get_p2() == pointlist[3] :
            l2 = line.Line(pointlist[2],pointlist[3])
            excluded.append(l2)
            if True or not np.allclose(pointlist[1],pointlist[2],rtol=.005):
                covered = line.Line(pointlist[1],pointlist[2])

        return res,excluded,covered


    def intersects_line_strings(self,line_strings):
        """
        Use the Shapely package to return the intersection 
        of this circle with the given line_string. 
        The shapely package approximates a circle using a
        64 point polygon. 
        """
        def splitLineString(lineString):
            res = []
            coords = lineString.coords
            p0 = coords[0]
            for k in range(1,len(lineString.coords)):
                p1 = coords[k]
                res.append(line.Line(p0,p1))
                p0 = p1
            return res
            
        def defragmentLineString(lineString):
            coords = lineString.coords
            res = []
            p0 = coords[0]
            res.append(p0)
            for k in range(1,len(lineString.coords)):
                p1 = coords[k]
                if not np.allclose(p0,p1,rtol=.005):
                    res.append(p1)
                p0 = p1
            if len(res) > 1:
               ls = LineString(res)
            else:
                ls = None
            return ls

        def defragment(multilineString):
            """
            Fix the roundoff error caused by shapely polygon 
            approximation.
            """
            retval = []
            for l in multilineString:
                ls = defragmentLineString(l)
                if ls is not None:
                    retval.append(ls)
            return MultiLineString(retval)  

        excluded = line_strings.difference(self.get_geometry())
        included = line_strings.difference(excluded)
        if type(excluded) == LineString:
            excluded = MultiLineString([excluded])

        cover = []
        if type(included) == LineString:
            dls = defragmentLineString(included)
            if dls is not None:
                cover = splitLineString(dls)
        else:
            for ls in included:
                dls = defragmentLineString(ls)
                if dls is not None:
                    cover = cover + splitLineString(ls)

        return defragment(excluded),cover
        

    def intersects_lines(self,line_set):
        """
        compute intersection of line set with circle and return a lists
        of segments that are excluded and included in the circles.
        
        Returns two lists :

        - excluded_lines : A possibly empty list consisting of segments outside the circle
        - included_lines: A possibly empty list consisting of segments inside the circle.
        
        
        """
        excluded_lines = []
        included_set = []
        for line in line_set:
            b,l,included = self.collides(line)
            excluded_lines = excluded_lines + l
            if b:
                included_set.append(included)
        return excluded_lines,included_set

    
    def compute_polar_slice_area(self,lines):
        """
        Given a set of lines enclosed in this circle which form a slice, compute the area that is enclosed
        between the lines and the circle using numerical integration. The lines are assumed to form 
        a polar function i.e. a radial line from the center of the circle can only intersect with a single
        line in the collection
        
        Parameters:

        lines- a set of connected lines INSIDE the circles. Lines are assumed not to double back
                    on themselves. i.e. no switchback patterns.  Lines are all on one side of the center. 
                    lines are assumed to be connected and form a POLAR function (i.e. a function of r,theta).
                    Note that this method will NOT work (throws exception if the lines do not form a polar 
                    function).
        Returns:
            - area included in the section between the lines and circumference of the circle 
                (the slice area). 

        Raises:
            - Exception if the lines do not form a polar function.
            
        """
        def pol2cart(theta, rho):
            x = rho * np.cos(theta)
            y = rho * np.sin(theta)
            return [x, y]

        translated_lines = []
        # Move the circle to the origin.
        for l in lines:
            try:
                p1 = list(np.subtract(l.get_p1(), self.center))
                p2 = list(np.subtract(l.get_p2(), self.center))
                newline = line.Line(p1,p2)
                translated_lines.append(newline)
            except:
                pdb.set_trace()
        
        inside_area = 0
        circle_area = self.area()
        nslices = 360
        dtheta = 2*np.pi/nslices
        # Go around the circle, computing intersections with the translated lines.
        for i in range(0,nslices):
            #ray is the ray inside the circle to compute intersection.
            theta = dtheta*i
            ray = line.Line([0,0],pol2cart(theta,self.get_radius()))
            intersection = None
            found = False
            # dist is the square of the distance from the center
            # to the line.
            dist = self.get_radius()*self.get_radius()
            for l in translated_lines:
                b,p = l.intersection(ray)
                if b :
                    if found :
                        raise Exception("Contour is not a polar function")
                    d = (p[0]**2 + p[1]**2)
                    if d  < dist:
                        found = True
                        intersection = p
                        dist = d

            if not found:
                # the ray did not intersect with a line.
                slice_area = circle_area / float(nslices)
            else:
                # compute the area of the triangular wedge:
                slice_area =  (np.pi*dist / float(nslices)) 
            inside_area = inside_area + slice_area
        excess_area = circle_area - inside_area
        return excess_area
        

    def __hash__(self):
        return hash(str(self.get_center()) + str(self.r))

    def __eq__(self,other):
        return self.get_center() == other.get_center() and self.r == other.r

    def __repr__(self):
        r = {}
        r['center'] = list(self.get_center())
        r['radius'] = self.r
        return json.dumps(r)

    def to_map(self):
        r = {}
        r['center'] = list(self.get_center())
        r['radius'] = self.r
        return r
