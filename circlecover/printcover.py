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

def printAntennaCircleCover(testName,testCircle,cover,coverage_file,points_to_cover):
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
    convex_hull = list(MultiPoint(points_to_cover).convex_hull.exterior.coords)
    result["convex_hull"] = list(convex_hull)
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
    result["testName"] = testName
    result["possible_centers"] = possible_centers
    result["cover_centers"] =  centers
    result["antenna_aperture"]= antenna_angle
    result["angles"] =  angles
    result["indexes"] = indices
    result["detection_coverage_file"] = coverage_file
    result["test_name"] = testName
    result["interference_contour"] = interference_contour
    result["algorithm"] = "AntennaCover"
    sea_excess_area,land_excess_area = excessarea.compute_excess_area_for_antenna_cover(indices, 
            angles, centers, coverage_file, possible_centers, interference_contour)
    result["sea_excess_area"] = sea_excess_area
    result["land_excess_area"] = land_excess_area
    sensor_centers = list(set([c[0] for c in cover]))
    result["sensor_loc"] = sensor_centers
    result["sensor_count"] = len(sensor_centers)

    output_file = testName + "AntennaCover." + str(antenna_angle) + ".json"
    f = open(output_file,"w")
    to_write = json.dumps(result,indent=4)
    f.write(to_write)
    f.close()



def printCover(line_endpoints,cover,centers,min_separation,covered_segments,testName, algorithm):
    """
    Store the test results in a json formatted file for later viewing..
    """
    if algorithm == FIXED_RADIUS:
        output_file = testName + '_F.json'
    elif algorithm == VAR_RADIUS:
        output_file = testName + '_V.json'
    else:
        output_file = testName + '_A.json'

    f = open(output_file,"w")
    result = {}
    result["line_endpoints"] = line_endpoints
    lines = []

    p0 = line_endpoints[0]
    for i in range(1,len(line_endpoints)):
        p1 = line_endpoints[i]
        l = Line(p0,p1)
        lines.append(l)
        p0 = p1

    earea,carea = excessarea.compute_excess_area(cover,lines)

    newcover = []
    for c in cover:
        newcover.append(c.to_map())
    result["cover"] = newcover
    result["centers"] =  centers
    result["min_separation"] = min_separation
    result["excess_area"] = earea
    result["cover_area"] = carea
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
    to_write = json.dumps(result)
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
    convex_hull = result["convex_hull"]
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

    for p in convex_hull:
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
    line_endpoints = result['interference_contour']
    interference_linestring = LineString(line_endpoints)
    circ = interference_linestring
    plot_coords(ax,interference_linestring,RED)
    plot_line(ax,interference_linestring,YELLOW)
    possible_centers = result['possible_centers']
    cover_centers = result['cover_centers']
    centers_linestring = LineString(cover_centers)
    possible_centers_linestring = LineString(possible_centers)
    circ = circ.union(centers_linestring)
    plot_coords(ax,centers_linestring,GREEN)
    plot_line(ax,possible_centers_linestring,BLUE)
    angles = result['angles']
    indexes = result['indexes']
    # load the antenna pattern file.
    coverage_filename = result['detection_coverage_file']
    detection_coverage = antennacover.read_detection_coverage(coverage_filename)
    for i in range(0,len(indexes)):
        polygon = detection_coverage[indexes[i]][1]
        rotated_cover = antennacover.rotate(polygon,angles[i])
        rotated_translated_cover = antennacover.translate(rotated_cover,cover_centers[i])
        circ = circ.union(rotated_translated_cover)
        p = PolygonPatch(rotated_translated_cover, fc=GRAY, ec=GRAY, alpha=0.5, zorder=2)
        ax.add_patch(p)

    xmin,ymin,xmax,ymax = circ.bounds

    #xmin = float(circ.bounds[0])
    #ymin = float(circ.bounds[1])
    #xmax = float(circ.bounds[2])
    #ymax = float(circ.bounds[3])

    delta = max(abs(xmax -xmin),abs(ymax-ymin))
    
    ax.set_xlim([xmin,xmin+delta])
    ax.set_ylim([ymin,ymin+delta])
    antenna_aperture = result['antenna_aperture']
    land_excess_area = result["land_excess_area"]
    sea_excess_area = result["sea_excess_area"]
    sensor_count = result["sensor_count"]
    title = "Algorithm = " + "Antenna_Cover; Antenna_aperture_angle = "  + str(antenna_aperture) +\
            "\nland_excess_area = " + str(land_excess_area) + " sea_excess_area = " + str(sea_excess_area) +\
            "\nsensor_count = " + str(sensor_count)
    plt.suptitle(title)
    plt.gcf().canvas.set_window_title(result["testName"] +  "_Antenna_" + str(antenna_aperture))
    if os.path.dirname(fileName) != '':
        mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(fileName))
    else:
        mpl.rcParams["savefig.directory"] = os.chdir("./")
    plt.show()




def show_results(fileName):
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
    
    line_endpoints = result['line_endpoints']
    interference_linestring = LineString(line_endpoints)
    plot_coords(ax,interference_linestring,RED)
    plot_line(ax,interference_linestring,YELLOW)
    centers = result['centers']
    centers_linestring = LineString(centers)
    plot_coords(ax,centers_linestring,GREEN)
    plot_line(ax,centers_linestring,BLUE)

    covr = result['cover']
    cover = []
    cover_centers = []
    for c in covr:
        cover.append(circle.Circle(center=c['center'],radius=c['radius']))
        cover_centers.append(c['center'])

    circ = cover[0].get_geometry()
    circ = circ.union(interference_linestring)
    circ = circ.union(centers_linestring)

    if len(cover_centers) > 1:
        cover_centers_linestring = LineString(cover_centers)
        plot_coords(ax,cover_centers_linestring,BLACK)
    else:
        plot_point(ax,cover_centers[0],BLACK)
        circ = circ.union(Point(cover_centers[0]))

    print( "Bounds = " + str(circ.bounds) )

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

    title = "Algorithm = " + result["algorithm"] + "\nsea_excess_area = " + str(format_e(result["excess_area"])) + "\ncover_area = " + str(format_e(result["cover_area"]))
    
    plt.suptitle(title)

    plt.gcf().canvas.set_window_title(result["testName"] +  "_" + result["algorithm"])

    if os.path.dirname(fileName) != '':
        mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(fileName))
    else:
        mpl.rcParams["savefig.directory"] = os.chdir("./")

    
    plt.show()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process arguments.')
    parser.add_argument("-f",help = "Result file")
    parser.add_argument("-a",help = "Result file")
    parser.add_argument("-t",help = "Result file")
    args = parser.parse_args()
    if args.f  is not None:
        fileName = args.f
        show_results(fileName)
    elif args.a is not None:
        fileName = args.a
        show_results_for_antenna_cover(fileName)
    elif args.t is not None:
        fileName = args.t
        show_results_for_antenna_circle_cover(fileName)

    
