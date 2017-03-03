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
import bisect
import operator
import excessarea
import circle as ccle
from line import Line
import logging
import circlecover
import copy
from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely import affinity
from collections import namedtuple


logger = logging.getLogger("circlecover")

def distance(point1,point2):
    """
    Returns the distance between point1 and point2.
    """
    return math.sqrt((point1[0] - point2[0])**2  + (point1[1] - point2[1])**2)


def pol2cart(r,theta):
    """
    Translate from polar to cartesian coordinates.
    """
    return (r*math.cos(float(theta)/180*math.pi), r*math.sin(float(theta)/180*math.pi))


def read_detection_coverage(fileName):
    """
    Read a bunch of antenna patterns and return an array of arrays containing the results.
    Returns a set tuples. The first element of the tuple is the max reach of the antenna
    The second element of the tuple is a polygon that gives the antenna shape.
    """
    f = open(fileName)
    lines = f.readlines()
    polarLocations = eval(lines[0])
    results = []
    for i in range(1,len(lines)):
        coverage = eval(lines[i])
        thiscover = [pol2cart(coverage[j],polarLocations[j]) for j in range(0,len(polarLocations))]
        max_distance_covered = np.max(coverage) 
        coverage_polygon = Polygon(thiscover)
        # Attach it as a tuple of < max_distance_covered, coverage_polygon  >
        cp  = namedtuple("cp",["max_coverage", "lobe"])
        results.append(cp(max_distance_covered,coverage_polygon))
    # Sort the points by max extent. This allows us to do a binary search
    # to find the detection coverage.
    results = sorted(results,key=lambda x : x[0])
    return results




def angle(center,r):
    """ return the angle between r1,r2 and the horizontal axis"""

    delta_x = abs(r[0] - center[0])
    delta_y = abs(r[1] - center[1])
    ang = math.atan2(delta_y,delta_x)

    if r[0] >= center[0] and r[1] < center[1]:
        # 4th. quad
        return 2*math.pi -ang
    elif r[0] >= center[0] and r[1] > center[1]:
        # 1st quad
        return ang
    elif r[0] <= center[0] and r[1] > center[1]:
        # second quad
        return math.pi - ang
    elif r[0] <= center[0] and r[1] < center[1]:
        # third quad
        return (math.pi + ang)
    else:
        return ang


def rotate(polygon, rotation_angle):
    """
    Rotate a polygon to a rotation_angle.
    """
    return affinity.rotate(polygon,rotation_angle,origin=(0,0),use_radians = True)


def translate(polygon, point):
    """
    Translate a polygon from the origin to a point.
    """
    return affinity.translate(polygon,xoff=point[0],yoff=point[1])


def find_antenna_overlay_for_points(points_to_cover, center, radius, detection_coverage, antenna_angle):
    """
    Find the overlay pattern for antennas so they cover a sector.
    """
    def find_cover(rotated_antenna_pattern,points_to_cover):
        """ Find the subset of points covered by a rotated_antenna_pattern from a set of points_to_cover """
        return [point for point in points_to_cover if rotated_antenna_pattern.contains(Point(point))]

    def find_antenna_overlay_greedy(rotated_antenna_patterns,center,hull_to_cover,cover_patterns):
        """
        Do a greedy overlay of the points included in a circle.
        """
        max_cover = []
        max_pattern = None
        # do a greedy search to find that pattern.
        # We find the lobe orientation that consumes
        # the maximum porion of the hull.
        for pattern in rotated_antenna_patterns:
            # find the portion of the hull covered by the lobe.
            cover = find_cover(pattern[2],hull_to_cover)
            # track the maximum.
            if len(cover) > len(max_cover):
                max_cover = cover
                max_pattern = pattern

        if len(max_cover) == 0:
            # No more points left on the hull.
            return cover_patterns
        # Remove the points from the hull.
        for point in max_cover:
            hull_to_cover.remove(point)

        # add our greedy lobe to the cover_patterns
        cover_patterns.append(max_pattern)
        # remove it from the set we are considering
        rotated_antenna_patterns.remove(max_pattern)
        # recursively solve problem for remaining points
        return find_antenna_overlay_greedy(rotated_antenna_patterns,center,hull_to_cover,cover_patterns)
        
        
    # Compute convex hull of points to cover
    mp = MultiPoint(points_to_cover).convex_hull
    if isinstance(mp,Polygon):
        convex_hull = list(mp.exterior.coords)
    else:
        convex_hull = list(mp.coords)

    furthest_point_tuple = max([(p,distance(center,p)) for p in convex_hull], key = lambda t:t[1])
    radius = furthest_point_tuple[1]
    min_angle = min([angle(center,p) for p in convex_hull])
    # note - second element of tuple is None. (Not a syntax error)
    index = bisect.bisect_left(detection_coverage,(radius*1.2,))
    if index >= len(detection_coverage):
        print("radius " + str(radius*1.2))
        print("max_detection_coverage " + str(detection_coverage[len(detection_coverage) -1][0]))
        raise Exception("Antenna Pattern could not be found")
    # The unrotated antenna pattern
    antenna_pattern = detection_coverage[index].lobe
    # The number of patterns we allow around the circle.
    npatterns = int(2*math.pi / antenna_angle) + 1
    # the increment of each angle
    delta_angle = 2*math.pi / npatterns
    # Generate a set of rotated patterns.
    # These are rotated in increments of delta
    antenna_patterns_rotated = [(index,counter*delta_angle,translate(rotate(antenna_pattern,counter*delta_angle),center)) for counter in range(0,npatterns)]
    # get the greedy cover
    lobe_patterns =   find_antenna_overlay_greedy(antenna_patterns_rotated,center,convex_hull,[])
    lobe_covers = []
    # now track the points covered by each lobe.
    for pattern in lobe_patterns:
        lobe_cover = []
        for point in points_to_cover:
            if pattern[2].contains(Point(point)):
                lobe_cover.append(point)
        for point in lobe_cover:
            # Remove the point from the collection
            # of points to cover
            points_to_cover.remove(point)
        lobe_covers.append(lobe_cover)

    # return a list of tuples 
    l =  zip(lobe_patterns,lobe_covers)
    return [tuple(list(s1) + [s2]) for s1,s2 in l]
                
            
            


def min_antenna_area_cover_greedy(possible_centers, interference_contour, antenna_cover_file, antenna_angle,  min_center_distance=0):
    """
    Greedy cover with variable sized discs with additional knot points added inside the
    interference contour that should be covered.

    Parameters:

        possible_centers : Set of locations (points) where sensors may be placed (typically on the shore line).
        interference_contour: A set of points denoting the boundary of the area that must be protected.
        antenna_cover_file: A file (generated offline using propagation modeling) 
                which indicates the area covered by each antenna with the axis of the antenna horizontal.
        antenna_angle : The aperture angle for the antenna.
        min_center_distance: The minimum distance between center locations that is permissible.

    Returns:

        A tuple (cover,points) where :

        cover: A vector  of antenna coverage patterns which covers the entire area enclosed by
                the interference_contour and centered at possible_centers.
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

                1.  Find the minimum circle cover.

                2.  Cover each circle with antennas such that the circles are completely covered.


    """

    def intersect_circle_with_points(circle,point_set):
        """ Find the subset of points within 1.2 radius of the circle and return the intersection and difference"""
        newcircle = ccle.Circle(circle.get_center(), circle.get_radius()*1.2)
        points_to_include = [p for p in point_set if newcircle.inside(p)]
        for p in points_to_include:
            point_set.remove(p)
        return points_to_include,point_set

    
        

    centers = copy.copy(possible_centers)
    antenna_cover_patterns = read_detection_coverage(antenna_cover_file)
    # Find the min circle cover.
    cover,covered_point_sets = circlecover.min_area_cover_greedy(centers,interference_contour,min_center_distance)
    
    antenna_coverage = []

    antenna_lobe = namedtuple("AntennaLobeExtended",["center","index","angle","lobe","covered_points"])
    for i in range(0,len(cover)):
        points_to_cover = covered_point_sets[i]
        radius = cover[i].get_radius()
        center = cover[i].get_center()
        coverage = find_antenna_overlay_for_points(points_to_cover, center, radius, antenna_cover_patterns,float(antenna_angle)/float(180) * math.pi)
        # Return a set of triples - (center,index,angle) where index is the index into the coverage map.
        # c[0] is the index c[1] is the angle and c[2] the points covered by that lobe.
        antenna_coverage = antenna_coverage +  [antenna_lobe(center,c[0],c[1],c[2],c[3])  for (k,c) in enumerate(coverage)]

    # we now have the antenna coverage in descending order of size. Now see if we can eliminate lobes.
    # by moving points into larger lobes. First arrange in ascending order.
    antenna_coverage.reverse()
    for i in range(0,len(antenna_coverage)):
        points_covered = antenna_coverage[i].covered_points
        for j in range(i+1,len(antenna_coverage)):
            lobe = antenna_coverage[j].lobe
            points_to_move = []
            for p in points_covered:
                if lobe.contains(Point(p)):
                    points_to_move.append(p)
            for p in points_to_move:
                antenna_coverage[j].covered_points.append(p)
                antenna_coverage[i].covered_points.remove(p)


    # we want to eliminate lobes that have a very small number of points included
    # to eliminate the noise. If an antenna lobe covers less than tolerance number
    # of points, then we can eliminate it. We pick it to be 1/10 percent of the 
    # number of grid points (arbitrarily -- should be passed in as a parameter).
    tolerance = .005*sum([len(t.covered_points) for t in antenna_coverage])

    print "tolerance ",tolerance

    # Now eliminate empty lobes.
    antenna_coverage.reverse()
    retval = []  
    for i in range(0,len(antenna_coverage)):
        if len(antenna_coverage[i].covered_points) > tolerance:
            print "covered_points ", len(antenna_coverage[i].covered_points)
            retval.append((antenna_coverage[i].center,antenna_coverage[i].index,antenna_coverage[i].angle))
            
    # Return a list of tuples that indicates the coverage.
    return retval
        



