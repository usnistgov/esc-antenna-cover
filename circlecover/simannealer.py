from simanneal import Annealer
from shapely.geometry import Polygon
from shapely.geometry import Point
import antennacover
import excessarea
import math
import pdb
import random
import copy

class SimAnneal(Annealer):

    def covers(self,cover_polygons):
        """ 

        boolean to determine whether the collection of points is covered by the given set of lobes. 
        If this returns True then a valid cover has been found. We allow a small amount of 
        slop in the cover (.005%) to give the algorithm a chance to do better.

        """
        not_covered_count = 0
        # For each point in our grid, check to see if the point is covered by at least one lobe
        for point in self.points_to_check:
            covered = False
            for i in range(0,len(cover_polygons)):
                if cover_polygons[i].contains(point):
                    covered = True
                    break
            if not covered:
                not_covered_count = not_covered_count + 1

        ratio = float(not_covered_count) / float(len(self.points_to_check))
        return ratio < .005


    def energy(self):
        """ 
        The energy function for the simulated annealing. 
        This is called back by the simanneal function. 
        """
        cover_polygons = self.cover_polygons
        union = cover_polygons[0]
        for i in range(1,len(cover_polygons)):
            union = union.union(cover_polygons[i])
        hull = union.convex_hull
        cost = hull.area
        return cost


    def move(self):
        """
        The move function called back by the simanneal parent class.
        """
        while True:
            # Pick a random angle
            index = random.randint(0,len(self.indexes) - 1)
            sign = -1 if random.randint(0,1) == 0 else 1
            # Peturb the solution by one degree.
            delta_angle = 1*math.pi/180*sign
            new_angle = delta_angle + self.state[index]
            new_angles = copy.copy(self.state)
            new_angles[index] = new_angle
            # Generate the cover lobes.
            cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,new_angles,self.centers,self.detection_coverage)
            if (self.covers(cover_polygons)):
                # if the cover is good then keep the new angles
                self.state = new_angles
                # keep the cover polygons.
                self.cover_polygons = cover_polygons
                break


    def get_result(self):     
        """ 
        Prune unecessary antenna lobes and retrieve the result.
        """
        cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,self.state,self.centers,self.detection_coverage)
        # Now remove a polygon at a time and see if the cover criterion is met.
        indexes_to_remove = []
        for i in range(0,len(cover_polygons)):
            newcover = [cover_polygons[k] for k in range(0,len(cover_polygons)) if k != i and k not in indexes_to_remove]
            if self.covers(newcover):
                indexes_to_remove.append(i)
            
        print ("indexes_to_remove " + str(indexes_to_remove))
        centers = [self.centers[i] for i in range(0,len(self.centers)) if i not in indexes_to_remove]
        indexes = [self.indexes[i] for i in range(0,len(self.indexes)) if  i not in indexes_to_remove]
        angles =  [self.state[i] for i in range(0,len(self.state)) if i not in indexes_to_remove]
        return zip(centers,indexes,angles)


    def __init__(self,interference_contour, possible_centers, detection_coverage_file,cover):
        self.interference_contour = interference_contour
        self.possible_centers = possible_centers
        self.centers = [c[0] for c in cover   ]
        self.indexes = [c[1] for c in cover   ]
        self.state =  [c[2] for c in cover  ]
        # load the detection coverage file.
        self.detection_coverage = antennacover.read_detection_coverage(detection_coverage_file)
        # the polygons representing the antenna shapes (rotated and translated)
        self.cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,self.state,self.centers,self.detection_coverage)
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
        
        self.steps = 1000

        
