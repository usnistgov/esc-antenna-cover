# Obsolete file -- pending deletion. Unfortunately did not work out too well.
from shapely.geometry import Polygon
from shapely.geometry import Point
import antennacover
import excessarea
import math
import pdb
import random
import copy
from scipy.optimize import basinhopping


class Annealer:

    def covers(self,angles):
        """ boolean to determine whether the collection of points is covered by the given set of lobes. 
            If this returns True then a valid cover has been found. """
        not_covered_count = 0
        
        for point in self.points_to_check:
            covered = False
            for i in range(0,len(self.cover_polygons)):
                if self.cover_polygons[i].contains(point):
                    covered = True
                    break
            if not covered:
                not_covered_count = not_covered_count + 1
                
        ratio = float(not_covered_count) / float(len(self.points_to_check))
        print "ratio = ",ratio
        if float(not_covered_count) / float(len(self.points_to_check)) < .005:
            return True
        else:
            return False


    def acceptTest(self,**kwargs):
        """ Return true if our solution is valid and false otherwise. """
        angles = kwargs["x_new"]
        if  self.covers(angles) :
            self.cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,angles,self.centers,self.detection_coverage)
            self.angles = angles
            return True
        else:
            old_angles = kwargs["x_old"]
            self.cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,old_angles,self.centers,self.detection_coverage)
            return False

            

    def cost_function1(self,angles):
        cover_polygons = self.cover_polygons
        union = cover_polygons[0]
        for i in range(1,len(cover_polygons)):
            union = union.union(cover_polygons[i])
        hull = union.convex_hull
        cost = hull.area
        print "cost = ",cost
        return cost
                
    def cost_function2(self,angles):
        intersection_area = 0
        cover_polygons = self.cover_polygons
        for i in range(0,len(cover_polygons)):
            for j in range(i+1,len(cover_polygons)):
                if cover_polygons[i].intersects(cover_polygons[j]):
                    intersection = cover_polygons[i].intersection(cover_polygons[j])
                    intersection_area = intersection_area + intersection.area
        cost =  math.INF if intersection_area == 0 else 1.0/intersection_area
        print "cost = ",cost
        return cost

    def takeStep(self,angles):
        index = random.randint(0,len(self.indexes) - 1)
        sign = -1 if random.randint(0,1) == 0 else 1
        delta_angle = math.pi/180*2*sign
        new_angle = delta_angle + self.angles[index]
        new_angles = copy.copy(self.angles)
        new_angles[index] = new_angle
        self.angles = new_angles
        self.cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,new_angles,self.centers,self.detection_coverage)
        return new_angles


    def anneal(self):
        basinhopping(self.cost_function1,self.angles,accept_test = self.acceptTest,take_step = self.takeStep)
        return zip(self.centers,self.indexes,self.angles)


    def __init__(self,interference_contour, possible_centers, detection_coverage_file,cover):
        self.interference_contour = interference_contour
        self.possible_centers = possible_centers
        self.centers = [c[0] for c in cover   ]
        self.indexes = [c[1] for c in cover   ]
        self.angles =  [c[2] for c in cover  ]
        # load the detection coverage file.
        self.detection_coverage = antennacover.read_detection_coverage(detection_coverage_file)
        # the polygons representing the antenna shapes (rotated and translated)
        self.cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,self.angles,self.centers,self.detection_coverage)
        # The bounding polygon consisting of the shore and the interference contour.
        self.bounding_polygon = excessarea.generate_bounding_polygon(possible_centers,interference_contour)
        # Take the union of these shapes
        union = self.cover_polygons[0]
        for i in range(1,len(self.cover_polygons)):
            union = union.union(self.cover_polygons[i])
        minx,miny,maxx,maxy = self.bounding_polygon.bounds
        ndivisions = 100
        deltax = float(maxx - minx) / float(ndivisions)
        deltay = float(maxy - miny) / float(ndivisions)
        self.points_to_check = [Point (minx+i*deltax, miny+j*deltay) 
                                for i in range(0,ndivisions) 
                                    for j in range(0,ndivisions) 
                                        if self.bounding_polygon.contains(Point (minx+i*deltax, miny+j*deltay))]
        

        
