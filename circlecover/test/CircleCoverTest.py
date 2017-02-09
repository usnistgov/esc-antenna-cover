import unittest
import sys
import os
sys.path.insert(0,os.path.abspath('../'))
import circle
import circlecover
import line
import pdb
import random
import math
import logging
import numpy as np
from line import Line
import printcover
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString

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

FIXED_RADIUS = "FixedRadiusPointCover"
VAR_RADIUS = "VariableRadiusLineCover"
AREA_COVER = "AreaCover"


def total_area(circle_collection):
    total_area = 0
    for c in circle_collection:
        total_area = total_area + c.area()
    return total_area




class CircleCoverTest(unittest.TestCase):

    def setUp(self):
        pass



    def testMinimumCircleSetCoverForLineSetGreedy(self):
        print("testMinimumAreaCircleSetCoverForLineSetGreedy")
        line_endpoints = [[20,80],[50,70],[60,60],[80,50],[60,40],[80,30],[100,30],[90,20]]
        centers = [(10,70),(30,60),(50,55),(50,30),(60,20)]
        testName =  "testMinimumCircleSetCoverForLineSetGreedy"
        circ,covered_segments = circlecover.min_line_cover_greedy(centers,line_endpoints)
        printcover.printCover(line_endpoints,circ,centers,0,covered_segments,testName, VAR_RADIUS)
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

        circ,covered_points = circlecover.min_point_cover_greedy_with_fixed_discs(centers,line_endpoints)
        printcover.printCover(line_endpoints,circ,centers,0,covered_segments,testName, FIXED_RADIUS)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)


    def testMinimumCircleSetCoverForLineSetGreedy2(self):
        print("testMinimumCircleSetCoverForLineSetGreedy2")
        line_endpoints = [[20,80],[50,70]]
        centers = [(10,70),(40,60)]
        savedCentrs = list(centers)
        circ,segments = circlecover.min_line_cover_greedy(centers,line_endpoints)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        testName = "testMinimumCircleSetCoverForLineSetGreedy2"
        printcover.printCover(line_endpoints,circ,savedCentrs,0,segments,testName, VAR_RADIUS)

    def testMinimumCircleSetCoverForLineSetGreedy3(self):
        print("testMinimumCircleSetCoverForLineSetGreedy3")
        line_endpoints = [[20,80],[50,70]]
        centers = [(10,70),(30,60)]
        circ,segments = circlecover.min_area_cover_greedy(centers,line_endpoints)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        testName = "testMinimumCircleSetCoverForLineSetGreedy3"
        printcover.printCover(line_endpoints,circ,centers,0,[],testName, AREA_COVER)

        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

    def testMinimumCircleSetCoverForLineSetGreedyRandom(self):
        print("testMinimumCircleSetCoverForLineSetGreedyRandom")
        line_endpoints = []
        line_segments = []
        random.seed(0)
        start = [100,100]
        end = [200,120]

        p1 = start[0]
        nsteps = 200
        for i in range(1,nsteps):
            p2 = random.randint(start[1],end[1])
            line_endpoints.append((p1,p2))
            p1 = start[0] + float(end[0] - start[0])/float(nsteps) * i

        lineString = LineString(line_endpoints)
        centers = []
        x_coords= [ random.randint(100,200) for i in range(1,150)]
        sorted_xcoords = sorted(x_coords)

        for p1 in sorted_xcoords:
            p2 = random.randint(90,100)
            centers.append((p1,p2))

        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        circ,included = circlecover.min_line_cover_greedy(centers,line_endpoints,min_center_distance = 20)
        testName = "LineCoverGreedyRandom"
        printcover.printCover(line_endpoints,circ,centers,20,included,testName, VAR_RADIUS)

       
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        u = Polygon()

        for c in circ:
            u  = u.union(c.get_geometry())

        self.assertTrue(u.buffer(1).contains(lineString))

        # Check if the segments in the cover are distinct subsets.
        for i in range(0,len(included)):
            for k in range(i+1,len(included) -1 ):
                for l in included[i]:
                    for m in included[k]:
                        self.assertFalse(l == m)

        circ,included = circlecover.min_point_cover_greedy_with_fixed_discs(centers,line_endpoints,min_center_distance=20)
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)
        printcover.printCover(line_endpoints,circ,centers,20,[],testName, FIXED_RADIUS)


    def testEscCoverVB(self):
        print("testEscCoverVB")
        esc_loc_x = [1771380,1769310,1769790,1768380,1767390,1764690,1762020,1759920,1753110,1741950,1752210,1757010,1761870,1768230,1772820,1777110,1781610,1786920,1793220]
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

        circ,included = circlecover.min_line_cover_greedy(centers,line_endpoints,min_center_distance = 60)
        testName = "VirginiaBeach"
        printcover.printCover(line_endpoints,circ,centers,60,included,testName,VAR_RADIUS)
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        circ,included = circlecover.min_point_cover_greedy_with_fixed_discs(centers,line_endpoints,min_center_distance = 60)
        printcover.printCover(line_endpoints,circ,centers,60,[],testName,FIXED_RADIUS)
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        circ,included = circlecover.min_area_cover_greedy(centers,line_endpoints,min_center_distance = 60)
        printcover.printCover(line_endpoints,circ,centers,60,[],testName,AREA_COVER)
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)


    def testEscCoverSF(self):
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
        testName = "SanFrancisco"
        circ,included = circlecover.min_line_cover_greedy(centers,line_endpoints,min_center_distance = 60)
        printcover.printCover(line_endpoints,circ,centers,60,included,testName,VAR_RADIUS)
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)
        circ,included = circlecover.min_point_cover_greedy_with_fixed_discs(centers,line_endpoints,min_center_distance = 60)
        printcover.printCover(line_endpoints,circ,centers,60,[],testName,FIXED_RADIUS)
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)
        circ,included = circlecover.min_area_cover_greedy(centers,line_endpoints,min_center_distance = 60)
        printcover.printCover(line_endpoints,circ,centers,60,[],testName,AREA_COVER)
        for point in line_endpoints:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)



    def testEstuary(self):
        """
        Test for a deep estuary.
        """
        interference_contour = [(20,55),(35,65),(40,60),(45,65),(50,55)]
        centers = [(20,46),(25,30),(30,20),(40,15),(50,30),(60,50)]
        circ,included = circlecover.min_area_cover_greedy(centers,interference_contour,min_center_distance=0)
        testName = "Estuary"
        printcover.printCover(interference_contour,circ,centers,0,[],testName,AREA_COVER)
        for point in interference_contour:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        circ,included = circlecover.min_point_cover_greedy_with_fixed_discs(centers,interference_contour,min_center_distance=0)
        printcover.printCover(interference_contour,circ,centers,0,[],testName,FIXED_RADIUS)
        for point in interference_contour:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)

        circ,included = circlecover.min_line_cover_greedy(centers,interference_contour,min_center_distance=0)
        printcover.printCover(interference_contour,circ,centers,0,[],testName,VAR_RADIUS)
        for point in interference_contour:
            flag = False
            for c in circ:
                if c.inside(point):
                    flag = True
                    break
            self.assertTrue(flag)
