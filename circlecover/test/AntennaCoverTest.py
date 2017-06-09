import unittest
import sys
import os
sys.path.insert(0,os.path.abspath('../'))
import circle
import line
import pdb
import random
import math
import copy
import logging
import numpy as np
from line import Line
from circle import Circle
import printcover
import antennacover
import annealer
import simannealer
import excessarea
import operator
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString

# set up logging to file - see previous section for more details
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='debug.txt',
                    filemode='w')
logger = logging.getLogger('antennacover')
# create file handler which logs even debug messages
sh = logging.StreamHandler()
fh = logging.FileHandler('debug.txt')
fh.setLevel(logging.DEBUG)
sh.setLevel(logging.DEBUG)
logger.addHandler(fh)
logger.addHandler(sh)

ANTENNA_COVER = "AntennaCover"


def total_area(circle_collection):
    total_area = 0
    for c in circle_collection:
        total_area = total_area + c.area()
    return total_area


def m_to_km(km):
    return [m/1000.0 for m in km]



class AntennaCoverTest(unittest.TestCase):

    def setUp(self):
        pass


    def testAntennaCoverForCircle(self):

        center = (10,10)
        radius = 100
        testCircle = circle.Circle(radius=radius,center=center)

        points_to_cover=  [ tuple(map(operator.add, (r*math.cos(theta*math.pi/float(180)),r*math.sin(theta*math.pi/float(180))), center)) for theta in range(270,360) for r in range(0,radius) ]
        points_to_cover=  points_to_cover +  [ tuple(map(operator.add, (r*math.cos(theta*math.pi/float(180)),r*math.sin(theta*math.pi/float(180))), center)) 
                            for theta in range(0,90) for r in range(0,radius) ]
        pointCount = len(points_to_cover)
        #60 degrees antenna angle.
        antenna_angle = math.pi/3
        coverage_file = "DetectionCoverage_60deg.txt"
        detection_coverage = antennacover.read_detection_coverage(coverage_file)
        coverage = antennacover.find_antenna_overlay_for_points(points_to_cover, center, radius, detection_coverage, antenna_angle)

        printcover._testPrintAntennaCircleCover("AntennaCircleCover",testCircle,coverage,coverage_file,60,points_to_cover)

        flag = True
        not_covered = 0
        for point in points_to_cover:
            p = Point(point)
            found = False
            for i in range(len(coverage)):
                if coverage[i][2].contains(p) :
                    found = True
                    break
            if not found:
                not_covered += 1

        print "Not Covered ", not_covered

        # We tolerate a small miss rate in the points.
        self.assertTrue(not_covered < .01 * pointCount)
        coverage = float(not_covered)/float(pointCount)
        
        print "not covered fraction = ", coverage



    def testAntennaCoverForCircle1(self):

        center = (20,10)
        radius = 120
        testCircle = circle.Circle(radius=radius,center=center)

        # Convert from polar to cartesian points in the sector of interest.
        points_to_cover = [ tuple(map(operator.add, (r*math.cos(theta*math.pi/float(180)),r*math.sin(theta*math.pi/float(180))), center)) for theta in range(95,195) for r in range(0,radius) ]
        
        pointCount = len(points_to_cover)

        copy_of_points_to_cover = copy.copy(points_to_cover)
        #60 degrees antenna angle.
        coverage_file = "detection-coverage/ITMDetectionCoverage_60deg.json"
        detection_coverage = antennacover.read_detection_coverage(coverage_file)
        coverage = antennacover.find_antenna_overlay_for_points(points_to_cover, center, radius, detection_coverage)
        printcover._testPrintAntennaCircleCover("AntennaCircleCover1",testCircle,coverage,coverage_file,60,points_to_cover)
        flag = True
        not_covered = 0
        for point in copy_of_points_to_cover:
            p = Point(point)
            found = False
            for i in range(len(coverage)):
                if coverage[i][1].contains(p) :
                    found = True
                    break
            if not found:
                not_covered += 1

        print "Not Covered ", not_covered

        not_covered_fraction =  not_covered / float(len(copy_of_points_to_cover))
        print "Not covered fraction ", not_covered_fraction
        # We tolerate a small miss rate in the points.
        self.assertTrue(not_covered < .01 * len(copy_of_points_to_cover))
        


    def testVB(self):
        # Convert all units to Km
        esc_loc_x =  [x/1000.0 for x in  [1771380,1769310,1769790,1768380,1767390,1764690,1762020,1759920,1753110,1741950,1752210,1757010,1761870,1768230,1772820,1777110,1781610,1786920,1793220]]
        esc_loc_y =  [y/1000.0 for y in [1827030,1817070,1806990,1797090,1787100,1776840,1767270,1756950,1746690,1735050,1727220,1717290,1707360,1697370,1687320,1677450,1667400,1657350,1647360]]
        ship_loc_x = [x/1000.0 for x in [1847012,1844913,1845660,1834150,1823280,1811715,1807512,1806671,1810710,1807769,1817910,1822503,1827218,1823623,1828432,1842183,1846928,1852378,1858591]]
        ship_loc_y = [y/1000.0 for y in [1843636,1833617,1823583,1811442,1799284,1787072,1777140,1767066,1759078,1749183,1741311,1731358,1721401,1709309,1699318,1691518,1681523,1671542,1661589]]

        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)

        line_endpoints = []
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            line_endpoints.append(p)

        testName = "VirginiaBeach"
        min_ctr_dist = 60
        poly = excessarea.generate_bounding_polygon(centers,line_endpoints)
        coverage_file = "detection-coverage/ITMDetectionCoverage_60deg.json"
        cover = antennacover.min_antenna_area_cover_greedy(centers,poly,coverage_file,min_center_distance = min_ctr_dist)
        printcover.printAntennaCover(testName, poly, centers, cover,coverage_file,min_ctr_dist)


    def testSF(self):
        esc_loc_x = m_to_km([-2300850,-2297160,-2284680,-2283390,-2284800,-2289540,-2287620,-2287740,-2287620,-2291760,-2289540,-2283720,
                                -2279730,-2254320,-2252430,-2253120,-2256900,-2273160,-2273970,-2273910])
        esc_loc_y = m_to_km([1986840,1977120,1966620,1957680,1947570,1937730,1926720,1917720,1907880,1897830,1887360,1876560,
                                1867620,1852470,1843620,1833720,1824660,1817640,1807710,1797930])
        ship_loc_x = m_to_km([-2414875,-2401190,-2405002,-2406089,-2407670,-2402495,-2402056,-2400759,-2390693,-2394883,-2393517,
                                -2394351,-2393072,-2396507,-2393492,-2394127,-2399780,-2396297,-2387368,-2387021])
        ship_loc_y = m_to_km([2019979,2007266,2001267,1992909,1982826,1970152,1959565,1950068,1937330,1927338,1917048,1908072,
                                1899767,1892623,1883306,1873335,1864772,1852304,1839577,1829662])

        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)
        line_endpoints = []
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            line_endpoints.append(p)

        testName = "SanFrancisco"
        min_ctr_dist = 60
        poly = excessarea.generate_bounding_polygon(centers,line_endpoints)
        cover = antennacover.min_antenna_area_cover_greedy(centers,poly,"detection-coverage/ITMDetectionCoverage_60deg.json",min_center_distance = min_ctr_dist)
        printcover.printAntennaCover(testName, poly, centers, cover,"60deg",min_ctr_dist)



    def testEstuary(self):
        """
        Test for a deep estuary.
        """
        interference_contour = [(20,55),(35,65),(40,60),(45,65),(50,55)]
        possible_centers = [(20,46),(25,30),(30,20),(40,15),(50,30),(60,50)]
        min_ctr_dist = 0
        poly = excessarea.generate_bounding_polygon(possible_centers,interference_contour)
        cover = antennacover.min_antenna_area_cover_greedy(possible_centers,poly,"detection-coverage/ITMDetectionCoverage_60deg.json",min_center_distance=min_ctr_dist)
        testName = "Estuary"
        printcover.printAntennaCover(testName, poly, possible_centers, cover,"detection-coverage/ITMDetectionCoverage_60deg.json",min_ctr_dist)


    def testEstuary1(self):
        """
        Test for a deep estuary.
        """
        interference_contour = [(20,55),(35,65),(40,60),(45,65),(50,55)]
        possible_centers = [(20,46),(25,30),(30,20),(40,15),(50,30),(60,50)]
        min_ctr_dist = 0
        poly = excessarea.generate_bounding_polygon(possible_centers,interference_contour)
        cover = antennacover.min_antenna_area_cover_greedy(possible_centers,poly,"detection-coverage/ITMDetectionCoverage_90deg.json",min_center_distance=min_ctr_dist)
        testName = "Estuary"
        printcover.printAntennaCover(testName, poly, possible_centers, cover,"detection-coverage/ITMDetectionCoverage_90deg.json",min_ctr_dist)


    def testEstuary2(self):
        """
        Test for a deep estuary.
        """
        interference_contour = [(20,55),(35,65),(40,60),(45,65),(50,55)]
        possible_centers = [(20,46),(25,30),(30,20),(40,15),(50,30),(60,50)]
        min_ctr_dist = 0
        poly = excessarea.generate_bounding_polygon(possible_centers,interference_contour)
        cover = antennacover.min_antenna_area_cover_greedy(possible_centers,poly,"detection-coverage/ITMDetectionCoverage_120deg.json",min_center_distance=min_ctr_dist)
        testName = "Estuary"
        printcover.printAntennaCover(testName, poly, possible_centers, cover,"detection-coverage/ITMDetectionCoverage_120deg.json",min_ctr_dist)

    def testEstuaryAnneal(self):
        """
        Test simulated annealing.
        """
        interference_contour = [(20,55),(35,65),(40,60),(45,65),(50,55)]
        possible_centers = [(20,46),(25,30),(30,20),(40,15),(50,30),(60,50)]
        antennacover.NDIVISIONS = 200
        min_ctr_dist = 0
        coverage_file = "detection-coverage/ITMDetectionCoverage_60deg.json"
        testName = "Estuary"
        poly = excessarea.generate_bounding_polygon(possible_centers,interference_contour)
        cover = antennacover.min_antenna_area_cover_greedy(possible_centers,poly,coverage_file,min_center_distance=min_ctr_dist)
        printcover.printAntennaCover(testName, poly, possible_centers, cover,coverage_file,min_ctr_dist)
        annealr = simannealer.SimAnneal(poly, coverage_file,cover)
        annealr.anneal()
        testName = "EstuaryAnneal"
        improved_cover = annealr.get_result()
        printcover.printAntennaCover(testName, poly, possible_centers, improved_cover,coverage_file,min_ctr_dist)

    def testSFAnneal (self):
        esc_loc_x = m_to_km([-2300850,-2297160,-2284680,-2283390,-2284800,-2289540,-2287620,-2287740,-2287620,-2291760,-2289540,-2283720,
                                -2279730,-2254320,-2252430,-2253120,-2256900,-2273160,-2273970,-2273910])
        esc_loc_y = m_to_km([1986840,1977120,1966620,1957680,1947570,1937730,1926720,1917720,1907880,1897830,1887360,1876560,
                                1867620,1852470,1843620,1833720,1824660,1817640,1807710,1797930])
        ship_loc_x = m_to_km([-2414875,-2401190,-2405002,-2406089,-2407670,-2402495,-2402056,-2400759,-2390693,-2394883,-2393517,
                                -2394351,-2393072,-2396507,-2393492,-2394127,-2399780,-2396297,-2387368,-2387021])
        ship_loc_y = m_to_km([2019979,2007266,2001267,1992909,1982826,1970152,1959565,1950068,1937330,1927338,1917048,1908072,
                                1899767,1892623,1883306,1873335,1864772,1852304,1839577,1829662])

        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)
        interference_contour = []
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            interference_contour.append(p)

        testName = "SanFrancisco"
        min_ctr_dist = 60
        coverage_file = "detection-coverage/ITMDetectionCoverage_60deg.json"
        poly = excessarea.generate_bounding_polygon(centers,line_endpoints)
        cover = antennacover.min_antenna_area_cover_greedy(centers,poly,coverage_file,min_center_distance = min_ctr_dist)
        printcover.printAntennaCover(testName, poly, centers, cover,coverage_file,min_ctr_dist)
        testName = "SanFranciscoAnneal"
        annealr = simannealer.SimAnneal(poly,coverage_file,cover)
        annealr.anneal()
        improved_cover = annealr.get_result()
        printcover.printAntennaCover(testName, poly, centers, improved_cover,coverage_file,min_ctr_dist)

    def testVBAnneal(self):
        # Convert all units to Km
        esc_loc_x =  [x/1000.0 for x in  [1771380,1769310,1769790,1768380,1767390,1764690,1762020,1759920,1753110,1741950,1752210,1757010,1761870,1768230,1772820,1777110,1781610,1786920,1793220]]
        esc_loc_y =  [y/1000.0 for y in [1827030,1817070,1806990,1797090,1787100,1776840,1767270,1756950,1746690,1735050,1727220,1717290,1707360,1697370,1687320,1677450,1667400,1657350,1647360]]
        ship_loc_x = [x/1000.0 for x in [1847012,1844913,1845660,1834150,1823280,1811715,1807512,1806671,1810710,1807769,1817910,1822503,1827218,1823623,1828432,1842183,1846928,1852378,1858591]]
        ship_loc_y = [y/1000.0 for y in [1843636,1833617,1823583,1811442,1799284,1787072,1777140,1767066,1759078,1749183,1741311,1731358,1721401,1709309,1699318,1691518,1681523,1671542,1661589]]

        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)

        interference_contour = []
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            interference_contour.append(p)

        testName = "VirginiaBeach"
        min_ctr_dist = 60
        coverage_file = "detection-coverage/ITMDetectionCoverage_60deg.json"
        poly = excessarea.generate_bounding_polygon(centers,interference_contour)
        cover = antennacover.min_antenna_area_cover_greedy(centers,poly,coverage_file,min_center_distance = min_ctr_dist)
        printcover.printAntennaCover(testName, poly, centers, cover,coverage_file,min_ctr_dist)
        testName = "VirginiaBeachAnneal"
        annealr = simannealer.SimAnneal(poly,coverage_file,cover)
        annealr.anneal()
        improved_cover = annealr.get_result()
        printcover.printAntennaCover(testName, poly, centers, improved_cover,coverage_file,min_ctr_dist)

    def testEastCoastAnneal(self):
        with open ("InterfContour_EastCoast.txt", "r") as myfile:
            data=myfile.readlines()
        esc_loc_x = [x/1000.0 for x in  eval(data[0])]
        esc_loc_y = [x/1000.0 for x in  eval(data[1])]
        ship_loc_x = [x/1000.0 for x in  eval(data[2])]
        ship_loc_y = [x/1000.0 for x in  eval(data[3])]
        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)

        interference_contour = []
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            interference_contour.append(p)

        testName = "EastCoast"
        min_ctr_dist = 0
        coverage_file = "DetectionCoverage_60deg.txt"
        poly = excessarea.generate_bounding_polygon(centers,interference_contour)
        cover = antennacover.min_antenna_area_cover_greedy(centers,poly,coverage_file,min_center_distance = min_ctr_dist)
        printcover.printAntennaCover(testName, poly, centers, cover,coverage_file,min_ctr_dist)
        testName = "EastCoastAnneal"
        annealr = simannealer.SimAnneal(poly,coverage_file,cover)
        annealr.anneal()
        improved_cover = annealr.get_result()
        printcover.printAntennaCover(testName, poly, centers, improved_cover,coverage_file,min_ctr_dist)
        
    def testEastCoastAnneal1(self):
        with open ("InterfContour_EastCoast.txt", "r") as myfile:
            data=myfile.readlines()
        esc_loc_x = [x/1000.0 for x in  eval(data[0])]
        esc_loc_y = [x/1000.0 for x in  eval(data[1])]
        ship_loc_x = [x/1000.0 for x in  eval(data[2])]
        ship_loc_y = [x/1000.0 for x in  eval(data[3])]
        centers = []
        for i in range(0,len(esc_loc_x)):
            center = (esc_loc_x[i],esc_loc_y[i])
            centers.append(center)

        interference_contour = []
        for i in range(0,len(ship_loc_x)):
            p = (ship_loc_x[i],ship_loc_y[i])
            interference_contour.append(p)

        testName = "EastCoast"
        min_ctr_dist = 0
        coverage_file = "detection-coverage/ITMDetectionCoverage_60deg.json"
        poly = excessarea.generate_bounding_polygon(centers,interference_contour)
        cover = antennacover.min_antenna_area_cover_anneal(centers,poly,coverage_file,min_center_distance = min_ctr_dist)
        printcover.printAntennaCover(testName, poly, centers, cover,coverage_file,min_ctr_dist)

    
    def testDPA(self):
        with open ("CandidateESCLoc_East.txt", "r") as myfile:
            data=myfile.readlines()
        
        esc_loc_x = [x/1000.0 for x in eval(data[0])]
        esc_loc_y = [y/1000.0 for y in eval(data[1])]
        
        centers = zip(esc_loc_x,esc_loc_y)
        
        with open ("DPA1.txt", "r") as myfile:
            data=myfile.readlines()


        dpa_loc_x = [x/1000.0 for x in eval(data[0])]
        dpa_loc_y = [y/1000.0 for y in eval(data[1])]
        dpa_loc = zip(dpa_loc_x,dpa_loc_y)
        
        poly = Polygon(dpa_loc)
        
        coverage_file = "detection-coverage/ITMDetectionCoverage_90deg.txt"  
        
        min_ctr_dist = 0

        cover = antennacover.min_antenna_area_cover_greedy(centers,poly,coverage_file,min_center_distance = min_ctr_dist)
        testName = "DPA1"
        printcover.printAntennaCover(testName,poly,centers,cover,coverage_file,min_ctr_dist)

        
