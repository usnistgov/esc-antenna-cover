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
from line import Line




FIXED_RADIUS = "FixedRadiusPointCover"
VAR_RADIUS = "VariableRadiusLineCover"
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

def printCoverForMatlab(line_endpoints,cover,centers,min_separation,covered_segments,testName, algorithm):
    """
    Print cover for viewing in matlab. (it has some nice click and show features).
    This is useful for debugging.
    """

    if algorithm == FIXED_RADIUS:
        output_file = testName + '_F.m'
    else:
        output_file = testName + '_V.m'

    p0 = line_endpoints[0]
    lines = []

    for i in range(1,len(line_endpoints)):
        p1 = line_endpoints[i]
        l = Line(p0,p1)
        lines.append(l)
        p0 = p1
        

    circ = cover[0].get_geometry()
    for i in range(1,len(cover)):
        circ = circ.union(cover[i].get_geometry())
    
    for l in lines:
        circ = circ.union(l)
    
    print "bounds = ", circ.bounds

    xmin = int(circ.bounds[0])
    ymin = int(circ.bounds[1])
    xmax = int(circ.bounds[2])
    ymax = int(circ.bounds[3])

    
        
    X = []
    Y = []
    r = []
    

    f = open(output_file,"w")

    f.write("close all\n")

    f.write("X = " + str(X) + ";\n")
    f.write("Y = " + str(Y) + ";\n")
    f.write("r = " + str(r) + ";\n")
    for li in lines:
        p1 = [li.get_p1()[0], li.get_p2()[0]]
        p2 = [li.get_p1()[1], li.get_p2()[1]]
        f.write("line(" + str(p1) + "," + str(p2)  + ");\n")

    centers = sorted(centers, key = lambda t:centers[1])

    #p1 = [centers[0][0] , centers[0][1]]
    #for i in range(1,len(centers)):
    #    p2 = [centers[i][0],centers[i][1]]
    #    f.write("line( "  + str(p1) + ", " + str(p2) + ")\n")
    #    p1 = p2

    lims = [xmin,xmax,ymin,ymax]
    f.write("axis(" + str(lims) + ")\n")
    X = []
    Y = []
    r = []
    minradius = .005*min([c.get_radius() for c in cover])
    if minradius < 1 :
        minradius = 1
    for c in centers:
        X.append(c[0])
        Y.append(c[1])
        r.append(minradius)

    f.write("X = " + str(X) + ";\n")
    f.write("Y = " + str(Y) + ";\n")
    f.write("r = " + str(r) + ";\n")
    f.write("centers = [X',Y'];\n")
    f.write("viscircles (centers,r','Color','b');\n")

    # Draw the cover
    colors = ['b','r','g','k']
    for i in range(0,len(cover)):
        center = cover[i].get_center()
        radius = cover[i].get_radius()
        color = colors[i%len(colors)]
        f.write("viscircles (" + str(list(center)) + "," + str(radius) + "," + "'Color'" + "," + "'" + str(color) + "');\n")
        if len(covered_segments)  != 0:
            for li in covered_segments[i]:
                p1 = [li.get_p1()[0], li.get_p2()[0]]
                p2 = [li.get_p1()[1], li.get_p2()[1]]
                f.write("line(" + str(p1) + "," + str(p2)   + ",'Color'" + "," + "'" + str(color) + "');\n")
        
    X = []
    Y = []
    r = []
    for c in cover:
        X.append(c.get_center()[0])
        Y.append(c.get_center()[1])
        r.append(c.get_radius())
    f.write("X = " + str(X) + ";\n")
    f.write("Y = " + str(Y) + ";\n")
    f.write("centers = [X',Y'];\n")
    f.write("r = " + str(r) + ";\n")
    #f.write("viscircles (centers,r','Color', 'r');\n")

    # Draw the centers of the  cover
    r = []
    for c in cover:
        r.append(minradius)
    f.write("r = " + str(r) + ";\n")
    f.write("viscircles (centers,r','Color','g');\n")
    logger.debug( "Computing Excess area: ")
    earea,carea = circlecover.compute_excess_area(cover,lines)
    f.write("total_area = " + str(total_area(cover)) + ";\n")
    logger.debug( "total_area = " + str(total_area(cover)))
    logger.debug ("excess_area = "+ str(earea))
    logger.debug( "cover_area = " + str(carea))

    f.write("excess_area = " +  str(earea) + ";\n")
    f.write("min_separation = " + str(min_separation) + ";\n")
    f.write("title({'ALGORITHM : " + algorithm  + "','minSeparation = " + str(min_separation)   + "','totalArea = " + str(total_area(cover)) + "','coverArea = " + str(carea) + "',' excessArea = " + str(earea) + "'});")
    f.close()


def printCover(line_endpoints,cover,centers,min_separation,covered_segments,testName, algorithm):
    if algorithm == FIXED_RADIUS:
        output_file = testName + '_F.txt'
    else:
        output_file = testName + '_V.txt'

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
    def total_area(circle_collection):
        total_area = 0
        for c in circle_collection:
            total_area = total_area + c.area()

    def plot_coords(ax, ob, point_color):
        x, y = ob.xy
        ax.plot(x, y, 'o', color=point_color, zorder=1)

    f = open(fileName)
    result = json.load(f)
    plt.figure(dpi=90)
    # get the current axes.
    ax = plt.gca() 
    
    line_endpoints = result['line_endpoints']
    linestring = LineString(line_endpoints)
    plot_coords(ax,linestring,RED)
    centers = result['centers']
    centers_linestring = LineString(centers)
    plot_coords(ax,centers_linestring,GREEN)

    covr = result['cover']
    cover = []
    cover_centers = []
    for c in covr:
        cover.append(circle.Circle(center=c['center'],radius=c['radius']))
        cover_centers.append(c['center'])

    cover_centers_linestring = LineString(cover_centers)

    plot_coords(ax,cover_centers_linestring,BLACK)
    
        

    circ = cover[0].get_geometry()
    for i in range(1,len(cover)):
        circ = circ.union(cover[i].get_geometry())
    
    circ = circ.union(linestring)

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

    
