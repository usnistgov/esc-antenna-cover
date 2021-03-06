
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

"""
Excess area computation and miscellaneous utility routines used all over the place.
"""

from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
import circle
import line
import antennacover
import copy
import pdb
import numpy as np
import argparse
import json
import collections
from collections import namedtuple

def generate_bounding_polygon(possible_centers,interference_contour):
    """ generate a bounding polygon consisting of possible_centers and interference_contour points. """
    centers = copy.copy(possible_centers)
    # Create a multipoint polygon coonsisting of the original contour 
    # and the possible centers
    points = [point for point in interference_contour]
    # The centers and the shore points are listed in the same sorted order
    centers.reverse()
    for point in centers:
        points.append(point)
    # The polygon encloses the interference contour as well as the shore.
    # This is the region that needs to be covered.
    mp = Polygon(points)
    return mp

def round2(val):
    return round(val,2)

def generate_antenna_cover_polygons(indexes,angles,centers,detection_coverage):
    """ Generate the antenna cover rotated set of polygons """
    polygons = []
    for i in range(0,len(indexes)):
        polygon = detection_coverage[indexes[i]].lobe
        angle = angles[i]
        center = centers[i]
        rotated_translated_polygon = antennacover.translate(antennacover.rotate(polygon,angle),center)
        polygons.append(rotated_translated_polygon)
    return polygons

def compute_excess_area_for_antenna_cover(indexes, angles, centers, detection_coverage_file,
            possible_centers, protected_region,units="km"):
    """
    Parameters:
        indexes - the selected contours from detection_coverage_file
        angles - the orientation of each antenna.
        centers : the centers of the sensors (ordered).
        interference_contour : The interference contour point set (ordered)
        detection_coverage_file: The detection coverage file.
    """


    def point_covered_by_lobe(point,antenna_cover_lobes):
        for i in range(0,len(antenna_cover_lobes)):
            if antenna_cover_lobes[i].contains(point):
                return True
        return False
            


    # load the detection coverage file.
    detection_coverage = antennacover.read_detection_coverage(detection_coverage_file,coverage_units=units)
    
    # the bounding polygon representing the if-contour and the possible centers. This bounds the region that
    # needs to be covered.

    # Line string representing the interference contour. Add first and last point 
    # of possible_centers to this. This is the interference contour that is on the sea.
    interference_contour_linestring = LineString(protected_region)

    protected_polygon = Polygon(protected_region)
    
    # Line string representing the possible sensor locations 
    if len(possible_centers) > 1:
        possible_centers_linestring = LineString(possible_centers)
    else:
        possible_centers_linestring = Point(possible_centers[0])

    # the polygons representing the antenna shapes (rotated and translated)
    antenna_cover_polygons = generate_antenna_cover_polygons(indexes,angles,centers,detection_coverage)
    # Take the union of these shapes
    cover_union = antenna_cover_polygons[0]
    for i in range(1,len(antenna_cover_polygons)):
        cover_union = cover_union.union(antenna_cover_polygons[i])
    union = cover_union.union(interference_contour_linestring)
    minx,miny,maxx,maxy = union.bounds

    # Generate a point set and classify.
    ndivs = antennacover.NDIVISIONS
    deltax,deltay = float(maxx-minx)/ndivs, float(maxy-miny)/ndivs
    area_per_grid_point = deltax*deltay
    
    excess_sea_coverage_count = 0
    excess_land_coverage_count = 0
    outage_count = 0
    polygon_area_count = 0
    for i in range(0,ndivs):
        for j in range(0,ndivs):
            p = Point(minx + i*deltax, miny + j*deltay)
            if cover_union.contains(p) and not protected_polygon.contains(p):
                # The point p is now either on sea or land.
                d1 = interference_contour_linestring.distance(p)
                d2 = possible_centers_linestring.distance(p)
                if d1 <  d2 :
                    excess_sea_coverage_count = excess_sea_coverage_count + 1
                else:
                    excess_land_coverage_count = excess_land_coverage_count + 1
            if protected_polygon.contains(p):
                polygon_area_count = polygon_area_count + 1
                if not cover_union.contains(p):
                    outage_count = outage_count + 1
            

    # BUGBUG -- this does not always work because of BUGS with the data set (self intersecting polygons)
    #outage = protected_polygon.difference(union).area

    ExcessArea = namedtuple("ExcessArea","excess_sea_coverage excess_land_coverage outage_area protected_region")

    return ExcessArea(round2(excess_sea_coverage_count*area_per_grid_point), round2(excess_land_coverage_count*area_per_grid_point), round2(outage_count*area_per_grid_point), 
                    round2(polygon_area_count*area_per_grid_point))





def compute_excess_area_for_circle_cover(cover,possible_centers,protected_region):
    """
    compute the excess area for circle cover.
    Parameters:
        cover : A set of circles representing the circle cover of the region.
        interference_contour : The interference contour.
        possible_centers : The land loctions where the sensors may be placed.
    """

    # Line string representing the interference contour. Add first and last point 
    # of possible_centers to this. This is the interference contour that is on the sea.
    # we add the two points to complete the polygon
    interference_contour_linestring = protected_region.extererior.coords
    
    # Line string representing the possible sensor locations 
    possible_centers_linestring = LineString(possible_centers)

    cover_union = cover[0].get_geometry()
    for i in range(1,len(cover)):
        cover_union = cover_union.union(cover[i].get_geometry())
    # take the union of 
    union = cover_union.union(bounding_polygon)
    minx,miny,maxx,maxy = union.bounds

    # Generate a point set and classify.
    ndivs = antennacover.NDIVISIONS
    deltax,deltay = float(maxx-minx)/ndivs, float(maxy-miny)/ndivs
    area_per_grid_point = deltax*deltay
    
    excess_sea_coverage_count = 0
    excess_land_coverage_count = 0
    for i in range(0,ndivs):
        for j in range(0,ndivs):
            p = Point(minx + i*deltax, miny + j*deltay)
            if cover_union.contains(p) and not bounding_polygon.contains(p):
                # The point is inside the boundary of the cover but NOT in the bounding polygon
                # that needs to be covered.
                # The point p is now either on sea or land. We need to classify it.
                d1 = interference_contour_linestring.distance(p)
                d2 = possible_centers_linestring.distance(p)
                if d1 <  d2 :
                    excess_sea_coverage_count = excess_sea_coverage_count + 1
                else:
                    excess_land_coverage_count = excess_land_coverage_count + 1

    
    ExcessArea = namedtuple("ExcessArea","excess_sea_coverage excess_land_coverage outage_area protected_region")
    outage = bounding_polygon.difference(union).area
    protected_region = bounding_polygon.area
    return ExcessArea(round2(excess_sea_coverage_count*area_per_grid_point), round2(excess_land_coverage_count*area_per_grid_point),outage,protected_region)
                        
    

        
def compute_excess_area(circles, line_segments, grid_divisions=200):
    """
    This is a TESTING ONLY function.
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
        # Note that this ONLY works if the interference contour
        # is entirely to ONE side of the centers (as will usually 
        # be the case
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


if __name__=="__main__":
    """
    compute excess area from an antenna cover command line invocation utility.
    """
    def compute_excess_area_for_antenna_cover_from_result_file(fileName):
        with open(fileName) as f:
            result = json.load(f)
            ic_x = result["ic_x"]
            ic_y = result["ic_y"]
            interference_contour = [(ic_x[i],ic_y[i]) for i in range(0,len(ic_x))]
            esc_loc_x = result["esc_loc_x"]
            esc_loc_y = result["esc_loc_y"]
            possible_centers = [(esc_loc_x[i],esc_loc_y[i]) for i in range(0,len(esc_loc_x))]
            antenna_loc_x = result["antenna_loc_x"]
            antenna_loc_y = result["antenna_loc_y"]
            centers = [(antenna_loc_x[i],antenna_loc_y[i]) for i in range(0,len(antenna_loc_x))]
            detection_coverage_file = result['detection_coverage_file']
            indexes = result["indexes"]
            angles = result["angles"]
            points = zip(ic_x,ic_y)
            bounding_polygon = Polygon(points)
            return compute_excess_area_for_antenna_cover(indexes,angles,centers,
                        detection_coverage_file,possible_centers,bounding_polygon)

        
    def compute_excess_area_for_circle_cover_from_file(fileName):
        with open(fileName) as f :
            result = json.load(f)
            ic_x = result["ic_x"]
            ic_y = result["ic_y"]
            interference_contour = [(ic_x[i],ic_y[i]) for i in range(0,len(ic_x))]
            esc_loc_x = result["esc_loc_x"]
            esc_loc_y = result["esc_loc_y"]
            possible_centers = [(esc_loc_x[i],esc_loc_y[i]) for i in range(0,len(esc_loc_x))]
            sensor_loc_x = result["sensor_loc_x"]
            sensor_loc_y = result["sensor_loc_y"]
            radii = result["sensor_detection_radius"]
            cover = [circle.Circle(center=(sensor_loc_x[i],sensor_loc_y[i]),radius=radii[i]) for i in range(0,len(sensor_loc_x))]
            return compute_excess_area_for_circle_cover(cover,possible_centers,interference_contour)

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", default=None, help = "Result file for Antenna Cover ")
    parser.add_argument("-d", default=None, help = "Result file for Antenna Cover (DPA) ")
    parser.add_argument("-c", default=None, help = "Result file for circle cover")
    args = parser.parse_args()

    if not args.a and not args.c  and not args.d:
        print "Need -a or -c flag. Try -h to see options"
        exit()

    json_result = {}
    if args.a:
        filename = args.a
        result = compute_excess_area_for_antenna_cover_from_result_file(filename)
        json_result["excess_sea_coverage"] = round2(result.excess_sea_coverage)
        json_result["excess_land_coverage"] = round2(result.excess_land_coverage)
        json_result["outage_area"] = round2(result.outage_area)
        json_result["protected_area"] = round2(result.protected_region)
        print(str(json.dumps(json_result)))
    else:
        filename = args.c
        excess_sea_area,excess_land_area = compute_excess_area_for_circle_cover_from_file(filename)
        json_result["excess_sea_coverage"] = round2(result.excess_sea_coverage)
        json_result["excess_land_coverage"] = round2(result.excess_land_coverage)
        json_result["outage_area"] = round2(result.outage_area)
        json_result["protected_area"] = round2(result.protected_region)
        print(str(json.dumps(result)))
        

 

