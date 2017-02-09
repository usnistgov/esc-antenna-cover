
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
from shapely import affinity

logger = logging.getLogger("circlecover")

def distance(point1,point2):
    """
    Returns the distance between point1 and point2.
    """
    return math.sqrt((point1[0] - point2[0])**2  + (point1[1] - point2[1])**2)


def pol2cart(r,theta):
    return (r*math.cos(theta/180*math.pi), r*math.sin(theta/180*math.pi))


def read_antenna_cover_area(fileName):
    """
    Read a bunch of antenna patterns and return an array of arrays containing the results.
    Returns a set tuples. The first element of the tuple is the max reach of the antenna
    The second element of the tuple is a polygon that gives the antenna shape.

    """
    f = open(fileName)
    lines = readlines(f)
    polarLocations = eval(lines[0])
    results = []
    for i in range(1,len(lines)):
        coverage = eval(lines[i])
        thiscover = [pol2cart(coverage[j],polarlocations[j] for j in range(0,len(polarlocations))]
        maxcover = math.max(coverage)
        p = Polygon(thiscover)
        results.append((maxcover,p))
    return results




def angle(r1,r2):
    """ return the angle between r1,r2 and the horizontal axis"""
    dx = r2[1] - r1[1]
    dy = r2[0] - r1[0]
    if dx == 0 :
        return math.pi/2
    else:
        slope = dy/dx
        return math.atan(slope)

def rotate(polygon, rotation_angle):
    return affinity.rotate(polygon,rotation_angle,(0,0),use_radians = True)
        
    


def min_antenna_area_cover_greedy(possible_centers, interference_contour, antenna_cover_file,  min_center_distance=0):
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
                    that is the smallest antenna coverage shape for that point.
                    Let the set of such such coverage shapes be coverage_set.

                2. Find the biggest coverage area  C in coverage_set.

                3. Remove all points in interference_set that is enclosed by C.

                4. If the center of C exists in cover, increase the size of the existing circle.
                    otherwise add C to cover.

                5. If interference_set is empty, terminate.

                6. Remove all centers from the possible_centers set which violate the min_distance constraint.

                7. Recursively run Find_cover for interference_set




    """
    def find_tightest_enclosing_shape_for_points(centers,points):
        """
        Find the max-min circle. i.e. for each point of the collection
        find the tightest circle that encloses the point.
        Of all these circles return the one with the maximum diameter.
        """
        dist = [(center,point,distance(point,center)) for center in centers for point in points]
        retval = []
        # For each point, find the tightest shape that encloses the point. The center has to be
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
               max_min_angle = angle(r[1],r[2])

        return max_min_center[0],max_min_center[2],max_min_angle

    def find_cover(antenna_shape,points):
        """
        Find the set of points covered by the antenna shape.

        Parameters :
            antenna_shape - shape to test.
            points - set of points for which to check the cover.
`
        Note: We can do this in O(log (n)) time if we use a quad tree.

        """
        return [p for p in points if antenna_shape.inside(p)]


    def min_antenna_cover_greedy_worker(centers, interference_set,coverage_patterns, cover,points):
        # Find the max_min radius tightest enclosing circle.
        max_min_center, max_min_radius, rotation = find_tightest_enclosing_shape_for_points(centers,interference_set)

        # Check if the circle center aleady exists in our cover.
        found = False
        for c in cover:
            if c.get_center() == max_min_center:
                found = True
                break

        if found:
            cover.remove(c)

        # Circle center does not exist in our cover.
        # So add it to our cover.
        shape  = None
        for c in coverage_patterns:
            if c[0] > max_min_radius:
               shape = c[1]
               break
        # Rotate the shape to the axis coordinate.
        rotated_shape = rotate(shape,rotation)


        cover.append(rotated_shape)

        # find the points enclosed by the largest max_cover circle.
        max_cover = find_cover(rotated,interference_set)
        points.append(max_cover)
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
    ndivs = 100 
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
    logger.debug("interference_set size " + str(len(interference_set)))
    logger.debug("interference_contour size " + str(len(interference_contour)))
    antenna_cover = read_antenna_cover_area(antenna_cover_file)
    cover,covered = min_antenna_cover_greedy_worker(centers, interference_set, antenna_cover, cover,points)
    return cover,covered


