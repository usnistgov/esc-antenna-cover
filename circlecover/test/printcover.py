import sys
sys.path.append('../')
from shapely.geometry import Point
from shapely.geometry import LineString
from descartes import PolygonPatch
from matplotlib.patches import Polygon
from matplotlib import pyplot as plt
import circle
import pdb
import argparse
import json
import circlecover
import os
from line import Line

"""
Handles printing output of the tests and visualization of the output.
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



def printCover(line_endpoints,cover,centers,min_separation,covered_segments,testName, algorithm):
    """
    Store the test results in a json formatted file for later viewing..
    """
    if algorithm == FIXED_RADIUS:
        output_file = testName + '_F.txt'
    elif algorithm == VAR_RADIUS:
        output_file = testName + '_V.txt'
    else:
        output_file = testName + '_A.txt'

    if not os.path.exists("test-results") :
        os.mkdir("test-results")

    f = open("test-results/" + output_file,"w")
    result = {}
    result["line_endpoints"] = line_endpoints
    lines = []

    p0 = line_endpoints[0]
    for i in range(1,len(line_endpoints)):
        p1 = line_endpoints[i]
        l = Line(p0,p1)
        lines.append(l)
        p0 = p1

    earea,carea = circlecover.compute_excess_area(cover,lines)

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
    


def show_results(fileName):
    """
    Read the json formatted file previously stored and display it.
    This is to visualize the test results.
    """
    def total_area(circle_collection):
        total_area = 0
        for c in circle_collection:
            total_area = total_area + c.area()

    def plot_coords(ax, ob, point_color):
        x, y = ob.xy
        ax.plot(x, y, 'o', color=point_color, zorder=1)

    def plot_point(ax,point,point_color):
        ax.plot(point[0],point[1],'o',color=point_color,zorder=1)

    def plot_line(ax, ob, line_color):
        x, y = ob.xy
        ax.plot(x, y, color=line_color, alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)


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

    print "Bounds = ", str(circ.bounds)

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

    title = "Algorithm = " + result["algorithm"] + " excess_area " + str(result["excess_area"]) + " cover_area " + str(result["cover_area"])
    
    plt.suptitle(title)

    plt.show()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process arguments.')
    parser.add_argument("-f",help = "Result file")
    args = parser.parse_args()
    fileName = args.f
    show_results(fileName)

    
