import circle as ccle
import line
import itertools
import sys
import pdb
import numpy as np
import traceback
import math
import random
from line import Line
import logging
import copy
from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.geometry import Polygon

logger = logging.getLogger("circlecover")

def compute_excess_area(circles, line_segments, grid_divisions=200):
    """
    compute the excess area - i.e. the area between the line collection
    and the perephry of the circle using numerical integration.
    
    Returns a tuple -- the area and the grid size used for computation.
    The grid size is the length of each grid square.
    """

    def isoutside(circle, segments,point):
        """
        return true if the point is outside the line segment (between the 
        lines and the circumference of the circle.
        """
        l = line.Line(circle.get_center(),point)
        intersection_count = 0
        for line_segment in segments:
            if  line_segment.intersects(l):
                intersection_count = intersection_count + 1  
        #odd number of intersections means point is in 
        #excess area region.
        return (intersection_count % 2) == 1

    def get_line_segments_for_circle(circle,line_segments):
        """
        Compute the intersection of the circle with a set of line segments.
        """
        retval = []
        for l in line_segments:
            b,excluded,included = circle.collides(l)
            if b:
                retval.append(included)
        return retval


    def generate_grid(circles, grid_divisions):
        """
        Generate a grid for numerical integration.
        Returns a list of points enclosed in the circles.
        """
        centers = []
        left_bottom = None
        right_top = None
        for c in circles:
            if left_bottom == None:
                left_bottom = [float(c.get_center()[0] - c.get_radius()), float(c.get_center()[1] - c.get_radius())]
            else:
                if c.get_center()[0] - c.get_radius() < left_bottom[0]:
                    left_bottom[0] = c.get_center()[0] - c.get_radius()
                if c.get_center()[1] - c.get_radius() < left_bottom[1]:
                    left_bottom[1] = c.get_center()[1] - c.get_radius()
                    
            if right_top == None:
                right_top = [float(c.get_center()[0] + c.get_radius()), float(c.get_center()[1] + c.get_radius())]
            else:
                if c.get_center()[0] + c.get_radius() > right_top[0]:
                    right_top[0] = c.get_center()[0] + c.get_radius()
                if c.get_center()[1] + c.get_radius() > right_top[1]:
                    right_top[1] = c.get_center()[1] + c.get_radius()

        # find the side of the grid. We insist that each dimension should be broken up into
        # atleast grid_intervals intervals (arbitrary).
        x_distance = right_top[0] - left_bottom[0]
        y_distance = right_top[1] - left_bottom[1]

        if y_distance > x_distance:
            grid_size = y_distance/float(grid_divisions)
        else:
            grid_size = x_distance/float(grid_divisions)


        x_divisions = int(x_distance/grid_size)
        y_divisions = int(y_distance/grid_size)

        grid = []
        for i in range (0,x_divisions):
            for j in range (0,y_divisions):   
                point = [left_bottom[0] + i*grid_size,left_bottom[1] + j*grid_size]
                grid.append(point)

        # Find the subset of the grid that is inside at least one circle.
        newgrid = []
        for point in grid:
            found = False
            for c in circles:
                #if we find a circle that includes the point, 
                # then incude in our integration side.
                if not found and c.inside(point) :
                    newgrid.append(point)
                    found = True
                    
        # newgrid is a list of points within the circle set.
        # we can integrate over this grid. grid_size is the dimension of one side
        # of the integration square.
        return newgrid,grid_size

    grid,grid_size = generate_grid(circles,grid_divisions)
    grid_area = len(grid)*grid_size*grid_size

    count = 0
    for c in circles:
        segments = get_line_segments_for_circle(c,line_segments)
        # figure out the points in the annulus. Note that we consider a point to be in the annulus if 
        # the center is in corresponding grid point is in the annulus
        filteredPoints = filter(lambda p: c.inside((p[0] + grid_size/2, p[1] + grid_size/2)) and isoutside(c,segments,(p[0] + grid_size/2, p[1] + grid_size/2)),grid)
        for p in filteredPoints:
            grid.remove(p)
        count = count + len(filteredPoints)
    area = count * grid_size * grid_size
    return area,grid_area
        




def distance(point1,point2):
    """
    Returns the distance between point1 and point2.
    """
    return math.sqrt((point1[0] - point2[0])**2  + (point1[1] - point2[1])**2)

def min_area_cover_greedy(possible_centers, interference_contour, min_center_distance=0):
    """
    Greedy cover with variable sized discs with additional knot points added inside the
    interference contour that should be covered.
    
    """
    def find_tightest_enclosing_circle_for_points(centers,points):
        """
        Find the max-min circle. i.e. for each point of the collection
        find the tightest circle that encloses the point.
        Of all these circles return the one with the maximum diameter.
        """
        dist = [(center,point,distance(point,center)) for center in centers for point in points]
        retval = []
        # For each point, find the tightest circle that encloses the point. The center has to be
        # placed on one of the possible centers.
        for point in points:
            min_dist = min([ d for d in dist if d[1] == point ] , key = lambda t:t[2])
            retval.append(min_dist)

        max_min = 0
        max_min_center = None
        for r in retval:
            if r[2] > max_min:
               max_min_center = r
               max_min = r[2]
        
        return max_min_center[0],max_min_center[2]

    def find_cover(circle,points):
        """
        Find the set of points covered by the circle.
        """
        return [p for p in points if circle.inside(p)]
                
            

    def min_area_cover_greedy_worker(centers, interference_set,cover,points):
        max_min_center, max_min_radius = find_tightest_enclosing_circle_for_points(centers,interference_set) 
        max_circle = ccle.Circle(center=max_min_center,radius=max_min_radius)
        cover.append(max_circle)
        max_cover = find_cover(max_circle,interference_set)
        points.append(max_cover)
        # Remove the points from our cover from the interference
        # contour. This leaves the others to be covered.
        for p in max_cover:
            interference_set.remove(p)
        # Remove the center that we used to construct the circle.
        try:
            centers.remove(max_circle.get_center())
        except:
            pdb.set_trace()
        # Of the remaining centers, remove those that violate our distance criterion.
        centers_to_remove = [k for k,cntr in enumerate(centers) if distance(cntr,max_circle.get_center()) < min_center_distance ]
        for k in sorted(centers_to_remove, reverse=True):
            del centers[k]
        # if we have covered every point in the interference set, then we are done.
        if len(interference_set) == 0:
            return cover,points
        else:
            return min_area_cover_greedy_worker(centers,interference_set,cover,points)
        

 
    centers = copy.copy(possible_centers)

    ifcontour = copy.copy(interference_contour)
    # Create a multipoint polygon coonsisting of the original contour 
    # and the possible centers
    points = []
    for point in interference_contour:
        points.append(point)

    centers.reverse()
    for point in centers:
        points.append(point)
    

    mp = Polygon(points)

    # mp now contains the points of the interference contour as well as the 
    # shore locations.

    xmin = float(mp.bounds[0])
    ymin = float(mp.bounds[1])
    xmax = float(mp.bounds[2])
    ymax = float(mp.bounds[3])

    interference_set = []

    for point in interference_contour:
        interference_set.append(point)

    # ndivs is the number of divsions to break up the range (xmin,ymin,xmax,ymax)
    ndivs = 20
    deltaX = (xmax - xmin)/ndivs
    deltaY = (ymax - ymin)/ndivs
    # We add additional points inside the contour to make sure that our area
    # is entirely covered.
    for i in range(1,ndivs):
        for j in range(1,ndivs):
            x = xmin + deltaX*i
            y = ymin + deltaY*j
            if mp.contains(Point(x,y)):
                interference_set.append((x,y))

    points = []
    cover = []
    print "interference_set ", len(interference_set)
    print "interference_contour ", len(interference_contour)
    cover,covered = min_area_cover_greedy_worker(centers, interference_set,cover,points)
    return cover,covered


def min_point_cover_greedy_with_fixed_discs(possible_centers, interference_contour, min_center_distance=0):
    """
    Greedy cover with fixed size discs.
    
    """
    def find_tightest_enclosing_circle_for_points(centers,points):
        """
        Find the max-min circle. i.e. for each point of the collection
        find the tightest circle that encloses the point.
        Of all these circles return the one with the maximum diameter.
        """
        dist = [(center,point,distance(point,center)) for center in centers for point in points]
        retval = []
        for point in points:
            min_dist = min([ d for d in dist if d[1] == point ] , key = lambda t:t[2])
            retval.append(min_dist)

        max_min = 0
        max_min_center = None
        for r in retval:
            if r[2] > max_min:
               max_min_center = r
               max_min = r[2]
        
        return max_min_center

    def find_center_with_max_cover(centers,points,radius):
        """
        Find the circle which covers the maximum number of points.
        """
        circles = [ccle.Circle(c,radius) for c in centers]
        max_circle = None
        max_cover = []
        for c in circles:
            cover = [p for p in points if c.inside(p)]
            if len(cover) > len(max_cover):
                max_circle = c
                max_cover = cover
        return max_circle, max_cover
                
            

    def min_cover_greedy_with_fixed_discs_worker(centers, radius,interference_contour,cover,points):
        max_circle,max_cover = find_center_with_max_cover(centers,interference_contour,radius)
        cover.append(max_circle)
        points.append(max_cover)

        # Remove the points from our cover from the interference
        # contour. This leaves the others to be covered.
        for p in max_cover:
            interference_contour.remove(p)

        centers.remove(max_circle.get_center())

        centers_to_remove = [k for k,cntr in enumerate(centers) if distance(cntr,max_circle.get_center()) < min_center_distance ]

        for k in sorted(centers_to_remove, reverse=True):
            del centers[k]

    

        if len(interference_contour) == 0:
            return cover,points
        else:
            if len(centers) == 0:
                print "Cover not found increase radius"
                return None,None
            return min_cover_greedy_with_fixed_discs_worker(centers,radius,interference_contour,cover,points)
        
    max_min_center, max_min_point, max_min_radius = find_tightest_enclosing_circle_for_points(possible_centers,interference_contour) 

    coverNotFound = True
    while coverNotFound:
        cover = []
        points = []
        centers = copy.copy(possible_centers)
        ifcontour = copy.copy(interference_contour)
        cover,covered =  min_cover_greedy_with_fixed_discs_worker(centers,max_min_radius,ifcontour,cover,points)
        if cover is None:
            max_min_radius = 1.1 * max_min_radius
            print "Retrying with radius = ", max_min_radius
        else:
            break
    
    return cover,covered
    


def min_line_cover_greedy(possible_centers, interference_contour, min_center_distance = 0):
    """
    Given a set of lines that form a multiline segment, cover the line segments
    with circles who's centers are closest to the intersection line.


    interference_contour : The interference contour.
    possible_centers : The possible centers where sensors may be placed.
    min_center_distance : The minimum distance between centers. Default value is 0.

    """

    def find_tightest_enclosing_circle(centers,cover,lines):
        max_min_distance =  -1e6
        max_min_center = None
        centers_to_check = list(set(centers + [c.get_center() for c in cover]))
        for line in lines:
            min_distance = 1e6
            mincenter = None
            # find the cheapest cover for each line. 
            for center in centers_to_check:
                # The circle radius that can cover both endpoints of the line.
                dist = line.circumscribing_radius(center)
                # find the smallest radius for this line
                if dist < min_distance:
                    min_distance = dist
                    mincenter = center
            # max_min_distance is the biggest minimum circle that covers any line in the
            # collection.
            if min_distance > max_min_distance:
                max_min_distance = min_distance
                max_min_center = mincenter

        return max_min_center,max_min_distance


    def min_line_cover_greedy_worker(centers,lines,cover,segments):

        # Compute the centers the closest centers

        max_min_center, max_min_radius = find_tightest_enclosing_circle(centers,cover,lines) 

        # Check if this circle is already in our cover.
            
        found = [(i,c) for (i,c) in enumerate ([ k for k in cover if max_min_center == k.get_center()])]
                
        if len(found) != 0:
            c = found[0][1]
            index = found[0][0]
            assert c.get_radius() <= max_min_radius
            newRadius = max(max_min_radius,c.get_radius())
            c.set_radius(newRadius)
            newlines,included_lines = c.intersects_lines(lines)
            segments[index] = segments[index] + included_lines
        else:
            c = ccle.Circle(max_min_center,max_min_radius)
            cover.append(c)
            newlines,included_lines = c.intersects_lines(lines)
            segments.append(included_lines)


        if len(newlines) == 0 :
            # no more line segments left to cover -- we are done.
            return cover,segments
        else:
            # Remove all the centers within our distance constraint so our constraint
            # is satisfied.
            if max_min_center in centers:
                 centers.remove(max_min_center)

        # eliminate all centers that are closer than the desired threshold
        centers_to_remove = [cntr for cntr in centers if cntr != max_min_center and distance(cntr,max_min_center) < min_center_distance ]
        for c in centers_to_remove:
            centers.remove(c)
        return min_line_cover_greedy_worker(centers,newlines,cover,segments)

    cover = []
    included_segments = []
    line_segments = []
    # Copy this because the list gets destroyed by the worker. We want the entries to be tuples.
    newcenters = [tuple(center) for center in possible_centers]
    # Construct a bunch of line segments corresponding to the interference contour
    p0 =  interference_contour[0]
    for i in range(1,len(interference_contour)):
        p1 = interference_contour[i]
        l = line.Line(p0, p1)
        line_segments.append(l)
        p0 = p1

    cloned_line_segments = copy.copy(line_segments)
     
    cover,covered = min_line_cover_greedy_worker(newcenters,line_segments,cover,included_segments) 

    
    return cover,covered
