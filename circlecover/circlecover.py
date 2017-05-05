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

import circle as ccle
import line
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



def min_area_cover_greedy(possible_centers, interference_contour, min_center_distance=0,ndivisions=100):
    """
    Greedy cover with variable sized discs with additional knot points added inside the
    interference contour that should be covered. This is an O(n^2) algorithm in ndivisions so the grid size
    has to be kept coarse.

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
                    possible_centers that is center of the the smallest circle that 
                    encloses that point.  Let the set of such circles be circle_set.

                2. Find the biggest diameter circle C in circle_set.

                3. Remove all points in interference_set that is enclosed by C.

                4. If the center of C exists in cover, increase the size of the existing circle.
                    otherwise add C to cover.

                5. If interference_set is empty, terminate.

                6. Remove all centers from the possible_centers set which violate the min_distance constraint.

                7. Recursively run Find_cover for interference_set




    """

    def generate_interference_set(possible_centers, interference_contour,ndivs):
        # Make a copy of the interference contour so we will not destroy our parameter list.
        ifcontour = copy.copy(interference_contour)
        # Create a multipoint polygon coonsisting of the original contour 
        # and the possible centers
        points = [point for point in interference_contour]
        # The centers and the shore points are listed in the same sorted order
        centers = copy.copy(possible_centers)
        centers.reverse()
        for point in centers:
            points.append(point)
        # The polygon encloses the interference contour as well as the shore.
        mp = Polygon(points)
        # We will fill the polygon with a grid and cover all the points in
        # the grid.
        xmin = float(mp.bounds[0])
        ymin = float(mp.bounds[1])
        xmax = float(mp.bounds[2])
        ymax = float(mp.bounds[3])

        interference_set = []

        for point in interference_contour:
            interference_set.append(point)
        # ndivs is the number of divsions to break up the range (xmin,ymin,xmax,ymax)
        ndivs = ndivisions 
        deltaX = (xmax - xmin)/ndivs
        deltaY = (ymax - ymin)/ndivs
        # We add additional points inside the polygon to make sure that our area
        # is entirely covered.
        for i in range(1,ndivs):
            for j in range(1,ndivs):
                x = xmin + deltaX*i
                y = ymin + deltaY*j
                if mp.contains(Point(x,y)):
                    interference_set.append((x,y))
        return mp, interference_set

    def covers(cover,protected_polygon):
        union = cover[0].get_geometry()
        for i in range(1,len(cover)):
            union = union.union(cover[i].get_geometry())
        uncovered = protected_polygon.difference(union)
        return uncovered.area/protected_polygon.area < .005

    def eliminate_redundant_circles(cover,cover_polygon):
        indexes_to_remove = []
        for i in range(0,len(cover)):
            testcover = [ cover[j] for j in range(0,len(cover)) if j != i and j not in indexes_to_remove]
            if covers(testcover,cover_polygon):
                indexes_to_remove.append(i)
        newcover = [cover[i] for i in range(0,len(cover)) if i not in indexes_to_remove]
        return newcover

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


    def min_area_cover_greedy_worker(centers, interference_set):
        # Find the max_min radius tightest enclosing circle.
        cover = []
        while len(interference_set) != 0:
            max_min_center, max_min_radius = find_tightest_enclosing_circle_for_points(centers,interference_set)

            # Check if the circle center aleady exists in our cover.
            en = [i for (i,c) in enumerate(cover) if c.get_center() == max_min_center] 
            # find the points enclosed by the largest max_cover circle.
            max_circle = ccle.Circle(center=max_min_center,radius=max_min_radius)
            max_cover = find_cover(max_circle,interference_set)
            found = len(en) > 0
            if found:
                # We alredy have a circle centered at the required location
                index = en[0]
                # Fetch the circle
                max_circle = cover[index]
                assert max_circle.get_radius() <= max_min_radius
                # Set a new raadius
                max_circle.set_radius(max_min_radius)
            else:
                # We dont have a circle at the center we found
                # Append the new circle to our cover
                cover.append(max_circle)

            # Remove the points from our cover from the interference
            # contour. This leaves the others to be covered.
            for p in max_cover:
                interference_set.remove(p)

            # Of the remaining centers, remove those that violate our distance criterion.
            centers_to_remove = [k for k,cntr in enumerate(centers)
                                if distance(cntr,max_circle.get_center()) < min_center_distance
                                    and max_circle.get_center() != cntr]
            for k in sorted(centers_to_remove, reverse=True):
                del centers[k]

        # if we have covered every point in the interference set, then we are done.
        return cover

    # This is the outer function

    # Make a copy of the centers
    centers = copy.copy(possible_centers)

    cover_polygon, interference_set = generate_interference_set(possible_centers,interference_contour,ndivisions)
    # Now call the worker function to do the hard work.
    cover = min_area_cover_greedy_worker(centers, interference_set)
    # Return the set of the circles and the subset of points each one covers 
    # as computed by the algorithm. 
    cover = eliminate_redundant_circles(cover,cover_polygon)

    print "nsensors ", len(cover)

    return cover


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
        max_min_distance =  -10e6
        max_min_center = None
        centers_to_check = list(set(centers + [c.get_center() for c in cover]))
        for line in lines:
            min_distance = 10e6
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
