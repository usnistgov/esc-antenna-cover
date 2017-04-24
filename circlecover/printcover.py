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

import sys
sys.path.append('../')
from shapely.geometry import Point
from shapely.geometry import MultiPoint
from shapely.geometry import LineString
from descartes import PolygonPatch
from matplotlib.patches import Polygon
from matplotlib import pyplot as plt
import circle
import pdb
import argparse
import json
import circlecover
import antennacover
import os
import excessarea
from line import Line
import matplotlib as mpl

"""

Handles printing output of the tests to a file and visualization of the output.

"""

FIXED_RADIUS = "FixedRadiusPointCover"
VAR_RADIUS = "VariableRadiusLineCover"
AREA_COVER = "AreaCover"
BLUE = '#6699cc'
GRAY = '#999999'
DARKGRAY = '#333333'
YELLOW = '#ffcc33'
GREEN = '#339933'
RED = '#ff3333'
BLACK = '#000000'

def total_area(circle_collection):
    total_area = 0
    for c in circle_collection:
        total_area = total_area + c.area()
    return total_area

def plot_coords(ax, ob, point_color):
    x, y = ob.xy
    ax.plot(x, y, 'o', color=point_color, zorder=1)

def plot_point(ax,point,point_color):
    ax.plot(point[0],point[1],'o',color=point_color,zorder=1)

def plot_line(ax, ob, line_color):
    x, y = ob.xy
    ax.plot(x, y, color=line_color, alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)



def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def _testPrintAntennaCircleCover(testName,testCircle,cover,coverage_file,points_to_cover):
    result = {}
    angles = []
    indices = []
    result['center'] = testCircle.get_center()
    for c in cover:
        indices.append(c[0])
        angles.append(c[1])
    result['indexes'] = indices
    result['angles'] = angles
    result["detection_coverage_file"] = coverage_file
    result["points_to_cover"] = list(points_to_cover)
    output_file = testName + ".json"
    f = open(output_file,"w")
    to_write = json.dumps(result,indent=4)
    f.write(to_write)
    f.close()


def printAntennaCover(testName, interference_contour, 
                      possible_centers, cover, coverage_file, 
                      antenna_angle, min_separation):
    result = {}
    centers = [c[0] for c in cover   ]
    indices = [c[1] for c in cover   ]
    angles =  [c[2] for c in cover  ]
    ic_x = [interference_contour[i][0] for i in range(0,len(interference_contour))]
    ic_y = [interference_contour[i][1] for i in range(0,len(interference_contour))]
    esc_loc_x = [possible_centers[i][0] for i in range(0,len(possible_centers))]
    esc_loc_y = [possible_centers[i][1] for i in range(0,len(possible_centers))]
    sensor_loc_x = [centers[i][0] for i in range(0,len(centers))]
    sensor_loc_y = [centers[i][1] for i in range(0,len(centers))]

    result["testName"] = testName
    #result["possible_centers"] = possible_centers
    result["esc_loc_x"] = esc_loc_x
    result["esc_loc_y"] = esc_loc_y

    result["cover_centers"] =  centers
    result["antenna_aperture"]= antenna_angle
    result["angles"] =  angles
    result["indexes"] = indices
    result["detection_coverage_file"] = coverage_file
    result["test_name"] = testName
    #result["interference_contour"] = interference_contour
    result["ic_x"] = ic_x
    result["ic_y"] = ic_y
    result["algorithm"] = "AntennaCover"
    sensor_centers = list(set([c[0] for c in cover]))
    #result["sensor_loc"] = sensor_centers
    result["sensor_loc_x"] = sensor_loc_x
    result["sensor_loc_y"] = sensor_loc_y
    result["sensor_count"] = len(sensor_loc_x)

    output_file = testName + "AntennaCover." + str(antenna_angle) + ".json"
    f = open(output_file,"w")
    to_write = json.dumps(result,indent=4,sort_keys=True)
    f.write(to_write)
    f.close()



def printCover(interference_contour,cover,centers,min_separation,covered_segments,testName, algorithm):
    """
    Store the circle cover test results in a json formatted file for later viewing.
    """
    if algorithm == FIXED_RADIUS:
        output_file = testName + '_F.json'
    elif algorithm == VAR_RADIUS:
        output_file = testName + '_V.json'
    else:
        output_file = testName + '_A.json'

    f = open(output_file,"w")
    ic_x = [interference_contour[i][0] for i in range(0,len(interference_contour))]
    ic_y = [interference_contour[i][1] for i in range(0,len(interference_contour))]
    esc_loc_x = [centers[i][0] for i in range(0,len(centers))]
    esc_loc_y = [centers[i][1] for i in range(0,len(centers))]
    lines = []

    p0 = interference_contour[0]
    for i in range(1,len(interference_contour)):
        p1 = interference_contour[i]
        l = Line(p0,p1)
        lines.append(l)
        p0 = p1


    centers_x = [c.get_center()[0] for c in cover]
    centers_y = [c.get_center()[1] for c in cover]
    radii = [c.get_radius() for c in cover]
    
    result = {}
    # Note the cover element is just for readability
    result["esc_loc_x"]  = esc_loc_x
    result["esc_loc_y"] = esc_loc_y
    result["ic_x"] = ic_x
    result["ic_y"] = ic_y
    result["sensor_loc_x"] = centers_x
    result["sensor_loc_y"] = centers_y
    result["sensor_detection_radius"] = radii
    result["min_separation"] = min_separation
    if covered_segments is not None:
        cseg = []
        for cover in covered_segments:
            lseg = []
            for k in cover:
                lseg.append(k.get_coordinates())
            cseg.append(lseg)
        result["covered_segments"] = cseg
    result["testName"] = testName
    result["algorithm"] = algorithm
    to_write = json.dumps(result,indent=4,sort_keys=True)
    f.write(to_write)
    f.close()
    


    
def show_results_for_antenna_circle_cover(fileName) :
    plt.figure(dpi=90)
    # get the current axes.
    ax = plt.gca()
    f = open(fileName)
    result = json.load(f)
    coverage_filename = result["detection_coverage_file"]
    indexes = result["indexes"]
    angles = result["angles"]
    center = result["center"]
    points_to_cover = result["points_to_cover"]
    detection_coverage = antennacover.read_detection_coverage(coverage_filename)
    circ = None
    for k in range(0,len(indexes)):
        polygon = detection_coverage[indexes[k]][1]
        rotated_cover = antennacover.rotate(polygon,angles[k])
        rotated_translated_cover = antennacover.translate(rotated_cover,center)
        if circ is None:
            circ = rotated_translated_cover
        else:
            circ = circ.union(rotated_translated_cover)
        p = PolygonPatch(rotated_translated_cover, fc=GRAY, ec=GRAY, alpha=0.5, zorder=2)
        ax.add_patch(p)

    for p in points_to_cover:
        plot_point(ax,p,RED)
    
    xmin = float(circ.bounds[0])
    ymin = float(circ.bounds[1])
    xmax = float(circ.bounds[2])
    ymax = float(circ.bounds[3])
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    
    if os.path.dirname(fileName) != '':
        mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(fileName))
    else:
        mpl.rcParams["savefig.directory"] = os.chdir("./")

    plt.show()

def show_results_for_antenna_cover(fileName):
    plt.figure(dpi=90)
    # get the current axes.
    ax = plt.gca()
    f = open(fileName)
    result = json.load(f)

    ic_x = result["ic_x"]
    ic_y = result["ic_y"]
    interference_contour = [(ic_x[i],ic_y[i]) for i in range(0,len(ic_x))]
    interference_linestring = LineString(interference_contour)

    plot_coords(ax,interference_linestring,RED)
    plot_line(ax,interference_linestring,YELLOW)

    esc_loc_x = result["esc_loc_x"]
    esc_loc_y = result["esc_loc_y"]
    possible_centers = [(esc_loc_x[i],esc_loc_y[i]) for i in range(0,len(esc_loc_x))]

    sensor_loc_x = result["sensor_loc_x"]
    sensor_loc_y = result["sensor_loc_y"]
    cover_centers = [(sensor_loc_x[i],sensor_loc_y[i]) for i in range(0,len(sensor_loc_x))]

    centers_linestring = LineString(cover_centers)
    possible_centers_linestring = LineString(possible_centers)

    plot_coords(ax,centers_linestring,GREEN)
    plot_line(ax,possible_centers_linestring,BLUE)
    angles = result['angles']
    indexes = result['indexes']

    # load the antenna pattern file.
    coverage_filename = result['detection_coverage_file']
    detection_coverage = antennacover.read_detection_coverage(coverage_filename)
    print "plotting polygons .... "

    # Compute the bounds of the polygon cover for drawing.
    circ = interference_linestring
    circ = circ.union(centers_linestring)
    for i in range(0,len(indexes)):
        polygon = detection_coverage[indexes[i]][1]
        rotated_cover = antennacover.rotate(polygon,angles[i])
        rotated_translated_cover = antennacover.translate(rotated_cover,cover_centers[i])
        circ = circ.union(rotated_translated_cover)
        p = PolygonPatch(rotated_translated_cover, fc=GRAY, ec=GRAY, alpha=0.5, zorder=2)
        ax.add_patch(p)
    xmin,ymin,xmax,ymax = circ.bounds
    delta = max(abs(xmax -xmin),abs(ymax-ymin))
    ax.set_xlim([xmin,xmin+delta])
    ax.set_ylim([ymin,ymin+delta])

    antenna_aperture = result['antenna_aperture']
    sensor_count = result["sensor_count"]
    print "computing excess area ... "
    sea_excess_area,land_excess_area = excessarea.compute_excess_area_for_antenna_cover(indexes, 
            angles, cover_centers, coverage_filename, possible_centers, interference_contour)
    title = "Algorithm = " + "Antenna_Cover; Antenna_aperture_angle = "  + str(antenna_aperture) +\
            "\nland_excess_area = " + str(land_excess_area) + " sea_excess_area = " + str(sea_excess_area) +\
            "\nsensor_count = " + str(sensor_count)
    plt.suptitle(title)
    plt.gcf().canvas.set_window_title(result["testName"] +  "_Antenna_" + str(antenna_aperture))
    if os.path.dirname(fileName) != '':
        mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(fileName))
    else:
        mpl.rcParams["savefig.directory"] = os.chdir("./")
    print "rendering ..."
    plt.show()




def show_results_for_circle_cover(fileName):
    """
    Read the json formatted file previously stored and display it.
    This is to visualize the test results.
    """
    def total_area(circle_collection):
        total_area = 0
        for c in circle_collection:
            total_area = total_area + c.area()



    f = open(fileName)
    result = json.load(f)
    plt.figure(dpi=90)
    # get the current axes.
    ax = plt.gca() 
    
    esc_loc_x = result["esc_loc_x"]
    esc_loc_y = result["esc_loc_y"]
    ic_x = result["ic_x"]
    ic_y = result["ic_y"]
    interference_contour = [(ic_x[i],ic_y[i]) for i in range(0,len(ic_x))]
    interference_linestring = LineString(interference_contour)
    plot_coords(ax,interference_linestring,RED)
    plot_line(ax,interference_linestring,YELLOW)
    sensor_loc_x = result["sensor_loc_x"]
    sensor_loc_y = result["sensor_loc_y"]
    possible_centers = [(esc_loc_x[i],esc_loc_y[i]) for i in range(0,len(esc_loc_x))]
    centers_linestring = LineString(possible_centers)
    plot_coords(ax,centers_linestring,GREEN)
    plot_line(ax,centers_linestring,BLUE)
    sensor_radii = result["sensor_detection_radius"]


    cover = [circle.Circle(center=(sensor_loc_x[i],sensor_loc_y[i]),radius=sensor_radii[i]) for i in range(0,len(sensor_loc_x))]
    cover_centers = [(esc_loc_x[i],esc_loc_y[i]) for i in range(0,len(esc_loc_x))]
    cover_union = cover[0].get_geometry()
    for i in range(1,len(cover)):
        cover_union = cover_union.union(cover[i].get_geometry())

    # Form a large geometry object so we can get the bounds of the picture
    circ = cover[0].get_geometry()
    circ = circ.union(interference_linestring)
    circ = circ.union(centers_linestring)

    if len(cover_centers) > 1:
        cover_centers_linestring = LineString(cover_centers)
        plot_coords(ax,cover_centers_linestring,BLACK)
    else:
        plot_point(ax,cover_centers[0],BLACK)
        circ = circ.union(Point(cover_centers[0]))

    for i in range(1,len(cover)):
        circ = circ.union(cover[i].get_geometry())
    

    xmin = float(circ.bounds[0])
    ymin = float(circ.bounds[1])
    xmax = float(circ.bounds[2])
    ymax = float(circ.bounds[3])
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    

    for ob in cover:
        p = PolygonPatch(ob.get_geometry(), fc=GRAY, ec=GRAY, alpha=0.5, zorder=2)
        ax.add_patch(p)

    print "computing excess area ... "
    sea_excess_area,land_excess_area = excessarea.compute_excess_area_for_circle_cover(cover, possible_centers, interference_contour)
    cover_area = cover_union.area
    title = "Algorithm = " + result["algorithm"] + "\nsea_excess_area = " + str(format_e(sea_excess_area)) +\
                            "\nland_excess_area = " + str(format_e(land_excess_area)) +\
                            "\ncover_area = " + str(format_e(cover_area))
    
    plt.suptitle(title)

    plt.gcf().canvas.set_window_title(result["testName"] +  "_" + result["algorithm"])

    if os.path.dirname(fileName) != '':
        mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(fileName))
    else:
        mpl.rcParams["savefig.directory"] = os.chdir("./")

    
    plt.show()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process arguments.')
    parser.add_argument("-c",help = "Result file for circle cover ")
    parser.add_argument("-a",help = "Result file for antenna cover algorithm")
    parser.add_argument("-t",help = "Result file for test antenna circle cover (testing only)")
    args = parser.parse_args()
    if args.c  is not None:
        fileName = args.c
        show_results_for_circle_cover(fileName)
    elif args.a is not None:
        fileName = args.a
        show_results_for_antenna_cover(fileName)
    elif args.t is not None:
        fileName = args.t
        show_results_for_antenna_circle_cover(fileName)

    
