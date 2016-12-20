import unittest
import sys
sys.path.append('../')
import circle
import circles
import line
import pdb
import random
import math
import logging
import numpy as np
from line import Line

from shapely.geometry import Point

# set up logging to file - see previous section for more details
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='debug.txt',
                    filemode='w')
logger = logging.getLogger('circlecover')
# create file handler which logs even debug messages
sh = logging.StreamHandler()
fh = logging.FileHandler('debug.txt')
fh.setLevel(logging.DEBUG)
sh.setLevel(logging.DEBUG)
logger.addHandler(fh)
logger.addHandler(sh)


def total_area(circle_collection):
    total_area = 0
    for c in circle_collection:
        total_area = total_area + c.area()
    return total_area

def slice_area(circle_collection, covered_lines):
    total_area = 0
    for i in range(0,len(circle_collection)):
        total_area = total_area + circle_collection[i].compute_slice(covered_lines[i])
    return total_area

def printCover(line_endpoints,cover,centers,min_separation,covered_segments,output_file):

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
    earea,carea = circles.compute_excess_area(cover,lines)
    f.write("total_area = " + str(total_area(cover)) + ";\n")
    logger.debug( "total_area = " + str(total_area(cover)))
    logger.debug ("excess_area = "+ str(earea))
    logger.debug( "cover_area = " + str(carea))

    f.write("excess_area = " +  str(earea) + ";\n")
    f.write("min_separation = " + str(min_separation) + ";\n")
    f.write("title({'GREEDY MAX MIN cover: minSeparation = " + str(min_separation)   + "','totalArea = " + str(total_area(cover)) + "','coverArea = " + str(carea) + "',' excessArea = " + str(earea) + "'});")
    f.close()
    

class CirclesCoverTest(unittest.TestCase):

    def setUp(self):
        pass



    def testMinimumAreaCoverForLine(self):
        l = line.Line([20,15],[37,37])
        c = [circle.Circle(center=[20,20], radius=10),circle.Circle(center=[35,25],radius=10),circle.Circle(center=[35,35],radius=5),circle.Circle(center=[20,50],radius=10), 
            circle.Circle(center=[30,30], radius=5)]
        (b,r) = circles.min_area_cover_for_line_brute_force(c,l)
        self.assertTrue(b)

    def testMinimumAreaCoverForLineGreedy(self):
        l = line.Line([20,15],[37,37])
        c = [circle.Circle(center=[20,20], radius=10),circle.Circle(center=[35,25],radius=10),circle.Circle(center=[35,35],radius=5),circle.Circle(center=[20,50],radius=10), 
            circle.Circle(center=[30,30], radius=5)]
        (b,r) = circles.min_area_cover_for_line_greedy(c,l)
        self.assertTrue(b)

    def testMinimumAreaCoverForLineSetBruteForce(self):
        lines = [line.Line([15,45],[38,44]),line.Line([38,44],[58,38]),line.Line([58,38],[66,49])]
        cir = [ circle.Circle(center=[20,35], radius=13), 
                    circle.Circle(center=[35,40], radius=13), 
                    circle.Circle(center=[45,55], radius=8),
                    circle.Circle(center=[65,43], radius = 9),
                    circle.Circle(center=[53,35], radius = 8) ,
                    circle.Circle(center=[53,35], radius=10) ]
        b,c = circles.min_area_cover_for_lines_brute_force(cir,lines)
        self.assertTrue(b)

        print "circles " , c

        soln_x = []
        soln_y = []
        soln_r = []
        for circ in c:
            soln_x.append(circ.get_center()[0])
            soln_y.append(circ.get_center()[1])
            soln_r.append(circ.get_radius())
        print "TODO : some assertion checking goes here"
        print "solnX = ",soln_x
        print "solnY = ",soln_y
        print "solnR = ",soln_r

    

    def testMinimumCircleSetCoverForLineSetGreedy4(self):
        line_endpoints = [[20,80],[50,70],[60,60],[80,50],[60,40],[80,30],[100,30],[90,20]]
        centers = [(10,70),(30,60),(50,55),(50,30),(60,20)]
        savedCentrs = list(centers)
        circ,covered_segments = circles.min_cover_greedy(centers,line_endpoints,min_center_distance=30)
        #circ,covered_segments = circles.min_excess_area_cover_greedy(centers,line_segments)
        self.assertTrue(len(circ) == len(covered_segments))
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        for p in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(p):
                    flag = True

            self.assertTrue(flag)

        # Check if the segments in the cover are distinct subsets.
        for i in range(0,len(covered_segments)):
            for k in range(i+1,len(covered_segments) -1 ):
                for l in covered_segments[i]:
                    for m in covered_segments[k]:
                        self.assertFalse(l == m)

        printCover(line_segments,circ,savedCentrs,0,covered_segments,"testMinimumCircleSetCoverForLineSetGreedy.m")

    def testMinimumCircleSetCoverForLineSetGreedy(self):
        line_endpoints = [[20,80],[50,70],[60,60],[80,50],[60,40],[80,30],[100,30],[90,20]]
        centers = [(10,70),(30,60),(50,55),(50,30),(60,20)]
        savedCentrs = list(centers)
        print "testMinimumCircleSetCoverForLineSetGreedy"
        circ,covered_segments = circles.min_cover_greedy(centers,line_endpoints)
        printCover(line_endpoints,circ,savedCentrs,0,covered_segments,"testMinimumCircleSetCoverForLineSetGreedy.m")
        #circ,covered_segments = circles.min_cover_greedy(centers,line_segments)
        #printCover(line_segments,circ,savedCentrs,0,covered_segments,"testMinimumCircleSetCoverForLineSetGreedyHull.m")
        self.assertTrue(len(circ) == len(covered_segments))
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
    
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        # Check if the segments in the cover are distinct subsets.
        for i in range(0,len(covered_segments)):
            for k in range(i+1,len(covered_segments) -1 ):
                for l in covered_segments[i]:
                    for m in covered_segments[k]:
                        self.assertFalse(l == m)


    def testMinimumCircleSetCoverForLineSetGreedy2(self):
        line_endpoints = [[20,80],[50,70]]
        centers = [[10,70],[40,60]]
        savedCentrs = list(centers)
        line_segments = []
        circ,segments = circles.min_cover_greedy(centers,line_endpoints)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.

        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)


        printCover(line_segments,circ,savedCentrs,0,segments,"testMinimumCircleSetCoverForLineSetGreedy2.m")

    def testMinimumCircleSetCoverForLineSetGreedy3(self):
        line_endpoints = [[20,80],[50,70]]
        centers = [[10,70],[30,60]]
        savedCentrs = list(centers)
        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        circ,segments = circles.min_area_cover_greedy(centers,line_endpoints)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        printCover(line_segments,circ,savedCentrs,0,segments,"testMinimumCircleSetCoverForLineSetGreedy3.m")

        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

    def testMinimumCircleSetCoverForLineSetGreedyRandom(self):
        line_endpoints = []
        centers = []
        line_segments = []
        random.seed(0)
        start = [100,100]
        end = [200,120]

        p1 = start[0]
        nsteps = 200
        for i in range(1,nsteps):
            p2 = random.randint(start[1],end[1])
            line_endpoints.append([p1,p2])
            p1 = start[0] + float(end[0] - start[0])/float(nsteps) * i

        for i in range(1,150):
            p1 = random.randint(100,200)
            p2 = random.randint(90,100)
            centers.append((p1,p2))
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        circ,included = circles.min_cover_greedy(centers,line_endpoints,min_center_distance = 20)
        printCover(line_endpoints,circ,centers,20,included,"testMinimumCircleSetCoverForLineSetGreedyRandom.m")

        #circ,included = circles.min_area_cover_greedy(centers,line_endpoints,min_center_distance=20)
        #printCover(line_endpoints,circ,centers,20,included,"testMinimumCircleSetCoverForLineSetGreedyRandom.m")

        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        # The returned fragments are enclosed in the 
        # corresponding circles. Lets check for that.
        for i in range(0,len(circ)):
            c = circ[i]
            lines = included[i]
            # check if each line segment is inside a circle.
            for k in lines:
                self.assertTrue(c.encloses(k))

        # Check if the segments in the cover are distinct subsets.
        for i in range(0,len(included)):
            for k in range(i+1,len(included) -1 ):
                for l in included[i]:
                    for m in included[k]:
                        self.assertFalse(l == m)


    def testEscCover(self):
        esc_loc_x = [1771380,1769310,1769790,1768380,176739,1764690,1762020,1759920,1753110,1741950,1752210,1757010,1761870,1768230,1772820,1777110,1781610,1786920,1793220]
        esc_loc_y = [1827030,1817070,1806990,1797090,1787100,1776840,1767270,1756950,1746690,1735050,1727220,1717290,1707360,1697370,1687320,1677450,1667400,1657350,1647360]
        ship_loc_x = [1847012,1844913,1845660,1834150,1823280,1811715,1807512,1806671,1810710,1807769,1817910,1822503,1827218,1823623,1828432,1842183,1846928,1852378,1858591]
        ship_loc_y = [1843636,1833617,1823583,1811442,1799284,1787072,1777140,1767066,1759078,1749183,1741311,1731358,1721401,1709309,1699318,1691518,1681523,1671542,1661589]

        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)

        line_endpoints = []
        
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            line_endpoints.append(p)

        circ,included = circles.min_cover_greedy(centers,line_endpoints,min_center_distance = 60)
        printCover(line_endpoints,circ,centers,60,included,"testEscCoverVaBeach.m")


    def testEscCover1(self):
        esc_loc_x = [-2300850,-2297160,-2284680,-2283390,-2284800,-2289540,-2287620,-2287740,-2287620,-2291760,-2289540,-2283720,-2279730,-2254320,-2252430,-2253120,-2256900,-2273160,-2273970,-2273910]
        esc_loc_y = [1986840,1977120,1966620,1957680,1947570,1937730,1926720,1917720,1907880,1897830,1887360,1876560,1867620,1852470,1843620,1833720,1824660,1817640,1807710,1797930]
        ship_loc_x =[-2414875,-2401190,-2405002,-2406089,-2407670,-2402495,-2402056,-2400759,-2390693,-2394883,-2393517,-2394351,-2393072,-2396507,-2393492,-2394127,-2399780,-2396297,-2387368,-2387021]
        ship_loc_y = [2019979,2007266,2001267,1992909,1982826,1970152,1959565,1950068,1937330,1927338,1917048,1908072,1899767,1892623,1883306,1873335,1864772,1852304,1839577,1829662]

        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)
        line_endpoints = []
        
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            line_endpoints.append(p)

        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))

        circ,included = circles.min_cover_greedy(centers,line_segments,min_center_distance = 60)
        printCover(line_endpoints,circ,centers,60,included,"testEscCoverSanFrancisco.m")

    def testEscCoverFixedDiscs(self):
        esc_loc_x = [1771380,1769310,1769790,1768380,176739,1764690,1762020,1759920,1753110,1741950,1752210,1757010,1761870,1768230,1772820,1777110,1781610,1786920,1793220]
        esc_loc_y = [1827030,1817070,1806990,1797090,1787100,1776840,1767270,1756950,1746690,1735050,1727220,1717290,1707360,1697370,1687320,1677450,1667400,1657350,1647360]
        ship_loc_x = [1847012,1844913,1845660,1834150,1823280,1811715,1807512,1806671,1810710,1807769,1817910,1822503,1827218,1823623,1828432,1842183,1846928,1852378,1858591]
        ship_loc_y = [1843636,1833617,1823583,1811442,1799284,1787072,1777140,1767066,1759078,1749183,1741311,1731358,1721401,1709309,1699318,1691518,1681523,1671542,1661589]

        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)

        line_endpoints = []
        
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            line_endpoints.append(p)

        circ,included = circles.min_point_cover_greedy_with_fixed_discs(centers,line_endpoints,min_center_distance = 60)
        printCover(line_endpoints,circ,centers,60,[],"testEscCoverVaBeachFixedDisk.m")
