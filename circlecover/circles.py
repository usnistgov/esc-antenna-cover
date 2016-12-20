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
from shapely.geometry import LineString
from shapely.geometry import MultiLineString

logger = logging.getLogger("circlecover")


def covers_line(circles,l):
    """
    We are given a set of circles and ask for the subset of the
    circles that covers the line l.
    
    return true if the set of circles covers the line.
  
    Parameters:
        - circles the set of circles to test.
        - cover - the subset of the circles that covers l
    
    Returns:

        - a pair of <true,array[circle]> representing the cover that the
        algorithm found or  <false,None> if no such cover exists.

    """

    def covers_line_worker(circles,l,cover):
        if len(circles) ==  0:
            return False, None
        # Get the topmost circle on the list
        top = circles[0]
        # check if it collides with the given line sebment.
        # lines contains the pieces caused by the collision.
        t,lines,contained =  top.collides(l)
        # Construct a new list of circles by removing the topmost element of the list.
        newc = list(circles)
        newc.remove(top)
       
        #the circle did not collide with the segment.
        if not t:
            # Test if the rest of the circles cover the line segment
            return covers_line_worker(newc,l,cover)

        # If the circle did collide with the line segment, include it in the cover.
        cover.append(top)
        # Now check if there are remaining line segments. If there are none, the line
        # is completely enclosed .
        if len(lines) == 0:
            # completely enclosed.
            return True,cover
        if len(lines) == 1:
            # If there is one segment, get the cover for that segment
            remainder = lines[0]
            return covers_line_worker(newc,remainder,cover)
        elif len(lines) == 2:
            # If there are two one segments, get the cover for those segments.
            left = lines[0]
            right = lines[1]
            b1, newCover =  covers_line_worker(newc,left,cover) 
            if not b1:
                return False,None
            else:
                return covers_line_worker(newc,right,newCover)
    return covers_line_worker(circles,l,[])



def min_area_cover_for_line_brute_force(circles,l):
    """
    The minimum area cover of a line selected from a set of circles.
    This is a brute force algorithm. It just tries to cover
    the lines using a set of circles by permuting the order of covering
    using different orders.

    Here we are given a set of circles and want to find the minimum cost cover
    for a line using brute force search.
    """
    def cost(circle_list):
        """
        return the cost of a circle (its area)
        """
        cst = 0
        for c in circle_list:
            cst = cst +  c.r**2
        return cst

    permuted_circles = itertools.permutations(circles)
    currentMinCover = None
    currentMinCost =  sys.maxint
    found = False
    for t in permuted_circles:
        found,cover =  covers_line(t,l)
        if found and cost(cover) < currentMinCost:
           currentMinCover = cover
           currentMinCost = cost(cover)
    return found,currentMinCover

def min_area_cover_for_line_greedy(circles,l):
    """
    The minimum area cover of a line selected from a set of circles.

    Sort the circles by radius and cover with the smallest radius first.
    """
    circles.sort(key=lambda x : x.get_radius)
    return covers_line(circles,l)

def min_area_cover_for_lines_brute_force(circles,lines):
    """
    The minimum area cover for a set of lines given
    a set of circles.
    """
    # coverset is the set of circles that covers the
    # given line set.
    coverset = []
    for line in lines:
        found,cover = min_area_cover_for_line_brute_force(circles,line)
        if not found:
            return False,None
        else:
            for c in cover:
                for r in coverset:
                    if r.get_center() == c.get_center():
                        if r.get_radius() < c.get_radius():
                            # Pick the larger circle and add it to the
                            # cover set.
                            if c in cover:
                                cover.remove(c)
                            coverset.remove(r)
                            coverset.append(c)
            coverset = coverset + cover
    return True, list(set(coverset))

        
def min_area_cover_for_lines_greedy(circles,lines):
    """
    The minimum area cover for a set of lines given
    a set of circles.
    The approach is to just cover each line in the lines set
    using the given circles starting from the smallest one.
    Then, take the uninon of these covers.

    This is for the situation where you know the circles and their
    radii apriori and you want to return the minimum area cover
    The algorithm proceeds by starting with the smallest circle and
    covering the lines, then recursively apply to the remaining lines.
    
    """
    # coverset is the set of circles that covers the
    # given line set.
    coverset = []
    for line in lines:
        found,cover = min_area_cover_for_line_greedy(circles,line)
        if not found:
            return False,None
        else:
            for c in cover:
                for r in coverset:
                    if r.get_center() == c.get_center():
                        if r.get_radius() < c.get_radius():
                            # Pick the larger circle and add it to the
                            # cover set.
                            if c in cover:
                                # Need to check here to ensure
                                # c is stil in cover - otherwise remove fails.
                                cover.remove(c)
                            coverset.remove(r)
                            coverset.append(c)
            coverset = coverset + cover
    return True, list(set(coverset))
        

def compute_excess_area(circles, line_segments, grid_divisions=200):
    """
    compute the excess area - i.e. the area between the line collection
    and the perephry of the circle using numerical integration.
    
    Returns the area and the grid size used for computation.
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
        



def circles_intersects_lines(circles,lines):
    linesToTest = lines
    for circle in circles:
        linesToTest,included = circle.intersects_lines(linesToTest)
    return linesToTest

def distance(point1,point2):
    return math.sqrt((point1[0] - point2[0])**2  + (point1[1] - point2[1])**2)


def min_point_cover_greedy_with_fixed_discs(possible_centers, interference_contour, min_center_distance=0):
    """
    Greedy cover with fixed discs.
    """

    def find_tightest_enclosing_circle_for_points(centers,points):
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

        for p in max_cover:
            interference_contour.remove(p)

        centers.remove(max_circle.get_center())

        centers_to_remove = [k for k,cntr in enumerate(centers) if distance(cntr,max_circle.get_center()) < min_center_distance ]

        for k in sorted(centers_to_remove, reverse=True):
            del centers[k]
    

        if len(interference_contour) == 0:
            return cover,points
        else:
            return min_cover_greedy_with_fixed_discs_worker(centers,radius,interference_contour,cover,points)
        
    max_min_center, max_min_point, max_min_radius = find_tightest_enclosing_circle_for_points(possible_centers,interference_contour) 
    cover = []
    points = []

    centers = copy.copy(possible_centers)
    ifcontour = copy.copy(interference_contour)

    cover,covered =  min_cover_greedy_with_fixed_discs_worker(centers,max_min_radius,ifcontour,cover,points)
    return cover,covered
    


def min_cover_greedy(possible_centers, interference_contour, min_center_distance = 0):
    """
    Given a set of lines that form a multiline segment, cover the line segments
    with circles who's centers are closest to the intersection line.

    algorithm is the algorithm used to locate the next center and radius

    interference_contour : The interference contour.
    possible_centers : The possible centers where sensors may be placed.
    min_center_distance : The minimum distance between centers.

    """

    def find_tightest_enclosing_circle(centers,lines):
        max_min_distance =  -1e6
        max_min_center = None

        for line in lines:
            min_distance = 1e6
            mincenter = None
            # find the cheapest cover for each line. 
            for center in centers:
                # The circle that can cover both endpoints of the line.
                dist = line.circumscribing_radius(center)
                if dist < min_distance:
                    min_distance = dist
                    mincenter = center
            # max_min_distance is the biggest minimum circle that covers any line in the
            # collection.
            if min_distance > max_min_distance:
                max_min_distance = min_distance
                max_min_center = mincenter

        return max_min_center,max_min_distance


    def min_cover_greedy_worker(centers,lines,cover,segments):

        # Compute the centers that are closest to the perpendicular bisector of each line.
        max_min_center, max_min_radius = find_tightest_enclosing_circle(centers,lines) 

        # Check if this circle is already in our cover.
    
        found = [(i,c) for (i,c) in enumerate ([ k for k in cover if max_min_center[0] == k.get_center()])]
        
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
            centers_to_remove = [k for k,cntr in enumerate(centers) if distance(cntr,c.get_center()) < min_center_distance ]
            for k in sorted(centers_to_remove, reverse=True):
                del centers[k]
            return min_cover_greedy_worker(centers,newlines,cover,segments)

    cover = []
    included_segments = []
    line_segments = []
    newcenters = copy.copy(possible_centers)
    p0 =  interference_contour[0]
    for i in range(1,len(interference_contour)):
        p1 = interference_contour[i]
        l = line.Line(p0, p1)
        line_segments.append(l)
        p0 = p1
    
    cover,covered = min_cover_greedy_worker(newcenters,line_segments,cover,included_segments) 

    return cover,covered

    

def improve1(cover,lines):
    """
    The following algorithm may result in an improvement for certain geometries.

    Given set of circes and a segment cover, check if cover can be improved
    by moving segments around between circles. We check the intersection areas
    between circles and ask if we can improve the solution by moving the segments
    in the intersection area to the smaller of the two circles.
    """
    circles = [ccle.Circle(center=c.get_center(), radius=c.get_radius()) for c in cover]
    covered = [ [line.Line(l.get_p1(),l.get_p2()) for l in li] for li in lines]
    # Check if we can find a better fit.
    # Move lines into smaller circles if possible.
    for i in range(0,len(circles)-1):
        for j in range(i+1,len(circles)):
            # a set of indices for lines that are enclosed in a smaller dia circle.
            try:
               inds = [k for (k,c) in enumerate(covered[i]) if circles[i].overlaps(circles[j]) and  circles[j].intersects(c) or circles[j].encloses(c)  ]
               if len(inds) != 0:
                   # Enhance the covered set for the j'th circle.
                   covered[j] = covered[j] + [covered[i][k] for k in inds]
                   # resize the j'th circle.
                   center = circles[j].get_center()
                   # Enhance the covered set for the j'th circle.
                   covered[j] = covered[j] + [covered[i][k] for k in inds]
                   # resize the j'th circle.
                   center = circles[j].get_center()
                   radius = np.max([l.circumscribing_radius(center) for l in covered[j]])
                   circles[j].set_radius(radius)
                   for index in sorted(inds, reverse=True):
                       del covered[i][index]
                   center = circles[i].get_center()
                   if len(covered[i]) == 0:
                       radius = 0
                   else:
                       radius = np.max([max(distance(center,l.get_p1()),distance(center,l.get_p2())) for l in covered[i]])
                   circles[i].set_radius(radius)
            except:
                pdb.set_trace()
    inds = [k for (k,c) in enumerate(circles) if c.get_radius() == 0]
    for i in sorted(inds, reverse=True):
        del covered[i]
        del circles[i]
    return circles,covered



def min_area_cover_greedy(possible_centers,interference_contour, min_center_distance = 0, nsegments=1):
    """
    Greedy minimum area cover.

    See discussion on 

    http://stackoverflow.com/questions/40748412/minimun-area-geometric-cover-for-a-set-of-line-segments

    Find the worst line, i.e. the line requiring the largest additional
    circle area for its best circle center option with the corresponding
    line segment entirely in the circle.

    Construct the corresponding circle.

    Remove fully covered line segments from the set and cut partially
    covered segments at the circle boundary.  

    Remove the center of the newly constructed circle from the set of possible centers.
    
    Iterate until no more line segments remain.

    Parameters:

    possible_centers : a list of points where the circles can be placed.
    line_segments: a list of Line which need to be optimally covered.
    min_center_distance: distance constraint for minimum distance between centers (defaults to 0).
    nsegments : # of segments to break up each line (defaullts to 10).
    

    Returns:

    return a tuple consisting of two arrays

    -  list of circles that completely covers the given lines.
    -  a collection of line segments for the lines enclosed by the corresponding 
       circles.


    """
        
    def distance(p1,p2):
        return math.sqrt((p1[0]  - p2[0])**2 + (p1[1] - p2[1])**2)


    def split_into_lines(line_segment):
        retval = []

    def splitLineString(lineString):
        res = []
        coords = lineString.coords
        p0 = coords[0]
        for k in range(1,len(lineString.coords)):
            p1 = coords[k]
            res.append(Line(p0,p1))
            p0 = p1
        return res

    def find_best_center(centers,line_segments):
        """
        find the min of the max of the min of each pairwise line segment circle cover for that line segment.
        """
        max_min_distance =  -1e6
        max_min_center = None
        max_min_cover_line = None

        for line_segment in line_segments:
            lines = splitLineString(line_segment)
            for line in lines:
                min_distance = 1e6
                mincenter = None
                min_point = None
                min_line = None
                # find the cheapest cover for each line. 
                for center in centers:
                    dist = line.circumscribing_radius(center)
                    if dist < min_distance:
                        min_distance = dist
                        mincenter = center
                        min_line = line
                # max_min_distance is the biggest minimum circle that covers any line in the
                # collection.
                if min_distance > max_min_distance:
                   max_min_distance = min_distance
                   max_min_center = mincenter
                   max_min_cover_line = min_line

        return max_min_center,max_min_distance,max_min_cover_line
    
    def min_area_cover_greedy_worker(centers,lines,cover,segments):
        center,radius,max_min_cover_line = find_best_center(centers,lines)
        # Check if this circle is already in our cover.
        found = [(i,c) for (i,c) in enumerate ([ k for k in cover if center == k.get_center()])]
        if len(found) != 0:
            c = found[0][1]
            index = found[0][0]
            assert c.get_radius() <= max_min_radius
            newRadius = max(max_min_radius,c.get_radius())
            c.set_radius(newRadius)
            newlines,included_lines = c.intersects_line_strings(lines)
            segments[index] = segments[index] + included_lines
        else:
            c = ccle.Circle(center,radius)
            cover.append(c)
            newlines,included_lines = c.intersects_line_strings(lines)
            segments.append(included_lines)



        if newlines.is_empty :
            # no more line segments left to cover -- we are done.
            return cover,segments
        else:
            # If there are remaining line segments, iterate on the portion 
            # that is left.
            centers.remove(center)
            # Remove all the centers within our distance constraint so our constraint
            # is satisfied.
            centers_to_remove = [k for k,cntr in enumerate(centers) if distance(cntr,c.get_center()) < min_center_distance ]
            for k in sorted(centers_to_remove, reverse=True):
                del centers[k]
            return min_area_cover_greedy_worker(centers,newlines,cover,segments)

    cover = [] 
    included_segments = []
    split_line_segments = []
    newcenters = copy.copy(possible_centers)
    lineString = MultiLineString([LineString(interference_contour)])
    cover,covered = min_area_cover_greedy_worker(newcenters,lineString,cover,included_segments) 
    
    return cover,covered
