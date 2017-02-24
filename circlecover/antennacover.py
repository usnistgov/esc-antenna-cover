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
        cp  = namedtuple("cp",["max_coverage", "coverage_area"])
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
        retval = []
        for point in points_to_cover:
            if rotated_antenna_pattern.contains(Point(point)):
                retval.append(point)
        return retval

    def find_antenna_overlay_greedy(rotated_antenna_patterns,center,hull_to_cover,cover_patterns):
        """
        Do a greedy overlay of the points included in a circle.
        """
        max_cover = []
        max_pattern = None
        for pattern in rotated_antenna_patterns:
            cover = find_cover(pattern[2],hull_to_cover)
            if len(cover) > len(max_cover):
                max_cover = cover
                max_pattern = pattern

        if len(max_cover) == 0:
            return cover_patterns
    
        for point in max_cover:
            hull_to_cover.remove(point)

        cover_patterns.append(max_pattern)
        rotated_antenna_patterns.remove(max_pattern)
        return find_antenna_overlay_greedy(rotated_antenna_patterns,center,hull_to_cover,cover_patterns)
        
        

    convex_hull = list(MultiPoint(points_to_cover).convex_hull.exterior.coords)
    furthest_point_tuple = max([(p,distance(center,p)) for p in convex_hull], key = lambda t:t[1])
    radius = furthest_point_tuple[1]
    min_angle = min([angle(center,p) for p in convex_hull])
    index = bisect.bisect_left(detection_coverage,(radius*1.2,))
    if index >= len(detection_coverage):
        raise Exception("Antenna Pattern could not be found")
    
    antenna_pattern = detection_coverage[index][1]
    npatterns = int(2*math.pi / antenna_angle) + 1
    delta_angle = 2*math.pi / npatterns
    antenna_patterns_rotated = [(index,counter*delta_angle,translate(rotate(antenna_pattern,counter*delta_angle),center)) for counter in range(0,npatterns)]

    return  find_antenna_overlay_greedy(antenna_patterns_rotated,center,convex_hull,[])
                
            
            
    
    
    

        
# def find_antenna_overlay_for_points(points_to_cover, center, radius, detection_coverage, antenna_angle):
#    """
#    Find the overlay pattern for antennas so they cover a sector.
#    """
#    convex_hull = MultiPoint(points_to_cover).convex_hull
#    pdb.set_trace()
#    angles = [angle(center,p) for p in convex_hull.exterior.coords]
#    radius = max([distance(center,p) for p in convex_hull.exterior.coords])
#    index = bisect.bisect_left(detection_coverage,(radius*1.2,))
#    antenna_pattern = detection_coverage[index][1]
#    min_angle = min(angles)
#    max_angle = max(angles)
#    pdb.set_trace()
#    delta = abs(max_angle - min_angle)
#    number_of_antennas =  2 + int(float(delta)/float(antenna_angle))
#    coverage = []
#    d_delta = float(delta)/float(number_of_antennas)
#    for i in range(0,number_of_antennas + 1 ):
#        azimuth_angle = min_angle + d_delta/2 + i*d_delta
#        rotated_coverage = rotate(antenna_pattern,azimuth_angle)
#        coverage.append((index,azimuth_angle,rotated_coverage))
#    return coverage
    


def min_antenna_area_cover_greedy(possible_centers, interference_contour, antenna_cover_file, antenna_angle,  min_center_distance=0):
    """
    Greedy cover with variable sized discs with additional knot points added inside the
    interference contour that should be covered.

    Parameters:

        possible_centers : Set of locations (points) where sensors may be placed (typically on the shore line).
        interference_contour: A set of points denoting the boundary of the area that must be protected.
        antenna_cover_file: A file (generated offline using propagation modeling) 
                which indicates the area covered by each antenna with the axis of the antenna horizontal.
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

    
        

    centers = copy.copy(possible_centers)
    antenna_cover_patterns = read_detection_coverage(antenna_cover_file)
    # Find the min circle cover.
    cover,covered_point_sets = circlecover.min_area_cover_greedy(centers,interference_contour,min_center_distance)
    
    # For each of the covered sets, find an antenna cover. 
    
    antenna_coverage = []

    cvg = namedtuple("cvg",["center","index","angle"])
    for i in range(0,len(covered_point_sets)):
        points_to_cover = covered_point_sets[i]
        radius = cover[i].get_radius()
        center = cover[i].get_center()
        coverage = find_antenna_overlay_for_points(points_to_cover, center, radius, antenna_cover_patterns,float(antenna_angle)/float(180) * math.pi)
        # Return a set of triples - (center,index,angle) where index is the index into the coverage map.
        antenna_coverage = antenna_coverage +  [cvg(center,c[0],c[1])  for (k,c) in enumerate(coverage)]


    return antenna_coverage
        


