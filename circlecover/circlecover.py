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

def distance(point1,point2):
    """
    Returns the distance between point1 and point2.
    """
    return math.sqrt((point1[0] - point2[0])**2  + (point1[1] - point2[1])**2)

def compute_excess_area(circles, line_segments, grid_divisions=200):
    """
    This is a support function that is used to evaluate a cover algorithm.
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
        k = 0
        for p in filteredPoints:
            grid.remove(p)
            k = k + 1
        count = count + k
    area = count * grid_size * grid_size
    return area,grid_area


def min_area_cover_greedy(possible_centers, interference_contour, min_center_distance=0):
    """
    Greedy cover with variable sized discs with additional knot points added inside the
    interference contour that should be covered.

    Parameters:

        possible_centers : Set of locations (points) where sensors may be placed (typically on the shore line).
        interference_contour: A set of points denoting the boundary of the area that must be protected.
        min_center_distance: The minimum distance between center locations that is permissible.

    Returns:

        A tuple (cover,points) where :

        cover: A vector  of Circles which covers the entire area enclosed by
                the interference_contour and possible_centers.
        points: An vector of vectors where each vector entry of the inner vector
                contains a set of Points covered
                by the corresponding Circle in cover.

    Algorithm:

        Prepare the area:

            1. Construct a polygon consisting of the possible_centers
                and interference_contour.
            2. Add additional points inside the constructed polygon so that the
                entire area inside he polygon is covered.

            Let the points constructed in the previous set
            (including the interference_contour but not including possible_centers)
            be denoted by interference_set

        Construct the cover for interference_set:

        Find Cover:
            Inputs:

                1. interference_set : Set of interference points that we want to cover
                            (from the last step).

                2. possible_centers: Set of possible centers.

                3. min_distance: Minimum distance permissible between circle centers.

            Let cover be the empty set. run the algorithm Find_cover

            Algorithm Find_cover:

                1. For each point in interference_set, find the center from possible_centers
                    possible_centers that is the smallest circle that encloses that point.
                    Let the set of such circles be circle_set.

                2. Find the biggest diameter circle C in circle_set.

                3. Remove all points in interference_set that is enclosed by C.

                4. If the center of C exists in cover, increase the size of the existing circle.
                    otherwise add C to cover.

                5. If interference_set is empty, terminate.

                6. Remove all centers from the possible_centers set which violate the min_distance constraint.

                7. Recursively run Find_cover for interference_set




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

        Parameters :

            circle - circle to test.
            points - set of points for which to check the cover.
`
        Note: We can do this in O(log (n)) time if we use a quad tree.

        """
        return [p for p in points if circle.inside(p)]


    def min_area_cover_greedy_worker(centers, interference_set,cover,points):
        # Find the max_min radius tightest enclosing circle.
        max_min_center, max_min_radius = find_tightest_enclosing_circle_for_points(centers,interference_set)

        # Check if the circle center aleady exists in our cover.
        found = False
        for c in cover:
            if c.get_center() == max_min_center:
                c.set_radius(max(max_min_radius,c.get_radius()))
                max_circle = c
                found = True
                break

        # Circle center does not exist in our cover.
        # So add it to our cover.
        if not found:
            max_circle = ccle.Circle(center=max_min_center,radius=max_min_radius)
            cover.append(max_circle)

        # find the points enclosed by the largest max_cover circle.
        max_cover = find_cover(max_circle,interference_set)
        points.append(max_cover)
        # Remove the points from our cover from the interference
        # contour. This leaves the others to be covered.
        for p in max_cover:
            interference_set.remove(p)
        # Remove the center that we used to construct the circle.
        centers.remove(max_circle.get_center())
        # Of the remaining centers, remove those that violate our distance criterion.
        centers_to_remove = [k for k,cntr in enumerate(centers)
                                if distance(cntr,max_circle.get_center()) < min_center_distance
                                    and max_circle.get_center() != cntr]
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
    ndivs = 60
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
    logger.debug("interference_set " + str(len(interference_set)))
    logger.debug("interference_contour " + str(len(interference_contour)))
    cover,covered = min_area_cover_greedy_worker(centers, interference_set,cover,points)
    return cover,covered


def min_point_cover_greedy_with_fixed_discs(possible_centers, interference_contour, min_center_distance=0):
    """
    Greedy cover with fixed size discs. Find a cover (circle centers)
    and fixed circle diameter for a set of points using circles of
    a fixed diameter where circles can be centered. Both the center
    locations of the circles AND the fixed diameter of the circles needs
    to be determined.

    Inputs:

    possible_centers : centers where circles may be placed
        (these points correspond to locations on the shore where sensors may be placed).
    interference_contour: A set of points that needs to be covered.
    min_center_distance: lower bound on the spacing between circles.

    Returns:

    A tuple (cover,covered) where:

    - cover: is a set of fixed diameter circles that completely cover interference contour points.
    - covered: A disjoint set of points covered by each circle as computed by the algorithm.

    Algorithm:

    1. Find the circle diameter:
        - For each point i in interference_contour:
            For each center j in possible_centers:
                find the distance d[ij] between i and j
        - Find the maximum distance D of all the d[ij] previously computed.

    2. Find the center C in possible_centers such that the circle of diameter D
        covers the greatest number of points in the interference_contour.

    3. Place the circle at C. Let Q be the set of points covered by the circle at C

    4. Remove Q from interference_contour. If interference_contour is empty, STOP.

    5. Remove C from the possible_centers.

    6. Remove all centers from possible_centers that are closer than min_center_distance to C.

    7. If interference_contour is not empty and possible_centers is empty:
            Increase D by some fraction (10%).
            restore the original original interference_contour

       Go to step 2.

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
        if max_circle is None:
            return None,None
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
                print("Cover not found increase radius")
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
            print( "Retrying with radius = " + str( max_min_radius ) )
        else:
            break

    return cover,covered


def min_line_cover_greedy(possible_centers, interference_contour, min_center_distance = 0):
    """
    Given a set of lines that form a multiline segment, cover the line segments
    with circles.


    Parameters:

    interference_contour : The interference contour.
    possible_centers : The possible centers where sensors may be placed.
    min_center_distance : The minimum distance between centers. Default value is 0.

    Algorithm:


     1. Find the worst line, i.e. the line requiring the largest additional
     circle area for its best circle center option with the corresponding
     line segment entirely in the circle.

     2. Construct / extend the corresponding circle.

     3. Remove fully covered line segments from the set and cut partially
     covered segments at the circle boundary. Generate a new set of segments.
     Iterate until no more line segments remain.


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
