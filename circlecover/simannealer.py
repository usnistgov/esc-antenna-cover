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
    def generate_test_points(self):
        min_x,min_y,max_x,max_y = self.bounding_polygon.bounds
        ndivs = 100
        delta_x = abs(float(max_x - min_x) / ndivs)
        delta_y = abs(float(max_y - min_y) / ndivs)
        self.points_to_test = [Point(min_x + i*delta_x,min_y + j*delta_y) \
                                    for i in range(0,100) \
                                        for j in range(0,100) \
                                            if self.bounding_polygon.contains(Point(min_x + i*delta_x,min_y + j*delta_y))]
                                            
        
    def compute_uncovered_ratio1(self,cover_polygons):
        count = 0
        # For each point in our grid, check to see if the point is covered by at least one lobe
        union = cover_polygons[0]
        for i in range(1,len(cover_polygons)):
             union = union.union(cover_polygons[i])
        for point in self.points_to_test:
            if not union.contains(point):
                count = count + 1

        ratio = float(count) / float(len(self.points_to_test))

        return ratio
        
    def compute_uncovered_ratio(self,cover_polygons):
        if len(cover_polygons) == 0:
           return 1
        if not self.flag:
           union = cover_polygons[0]
           for i in range(1,len(cover_polygons)):
               union = union.union(cover_polygons[i])
           return self.bounding_polygon.difference(union).area/self.bounding_polygon.area
        else:
           return self.compute_uncovered_ratio1(cover_polygons)
        

    def covers(self,cover_polygons):
        """ 

        boolean to determine whether the collection of points is covered by the given set of lobes. 
        If this returns True then a valid cover has been found. We allow a small amount of 
        slop in the cover (.005%) to give the algorithm a chance to do better.

        """

        ratio = self.compute_uncovered_ratio(cover_polygons)
        return  ratio <= self.max_ratio


    def energy(self):
        """ 
        The energy function for the simulated annealing. 
        This is called back by the simanneal function. 
        Our energy function is the area of the convex hull enclosed
        by the antenna lobes. 
        """
        cover_polygons = self.cover_polygons
        union = cover_polygons[0]
        for i in range(1,len(cover_polygons)):
            union = union.union(cover_polygons[i])
        return union.area


    def move(self):
        """
        The move function called back by the simanneal parent class.
        """
        while True:
            # pick an intersecting polygon index
            #pindex = random.randint(0, len(self.intersecting_polygons)-1)

            # pindex = self.index_count % len(self.intersecting_polygons)
            # index = self.intersecting_polygons[pindex]

            # index = self.index_count % len(self.cover_polygons)
            index = random.randint(0, len(self.cover_polygons)-1)
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
            self.index_count = self.index_count + 1


    def get_result(self):     
        """ 
        Prune unecessary antenna lobes and retrieve the result.
        """
        if len(self.centers) > 1 :
            cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,self.best_state,self.centers,self.detection_coverage)
            # Now remove a polygon at a time and see if the cover criterion is met.
            indexes_to_remove = []
            # loosen up our cover criterion a little
            for i in range(0,len(cover_polygons)):
                newcover = [cover_polygons[k] for k in range(0,len(cover_polygons)) if k != i and k not in indexes_to_remove]
                if self.covers(newcover):
                    indexes_to_remove.append(i)
            print ("indexes_to_remove " + str(indexes_to_remove))

        else:
            indexes_to_remove = []
            
        centers = [self.centers[i] for i in range(0,len(self.centers)) if i not in indexes_to_remove]
        indexes = [self.indexes[i] for i in range(0,len(self.indexes)) if  i not in indexes_to_remove]
        angles =  [self.best_state[i] for i in range(0,len(self.state)) if i not in indexes_to_remove]
        return zip(centers,indexes,angles)


    def __init__(self,bounding_polygon,detection_coverage_file,cover,steps = 0,tol=.005,coverage_units="km"):
        """
        Parameters:
            bounding_polygon : A polygon that defines the area to be covered.
            possible_centers : An ordered set of points on the coastline that follows the intereference contour with the same restrictions
                                as the interference contour. In general, because the coastline does not loop back on itself, 
                                these constraints can be naturally satisifed by just sampling the coastline.
            detection_coverage_file: A text file with the detection coverage of the antenna. 
                                See https://github.com/usnistgov/circle-cover/raw/master/figures/DetectionCoverage_90deg.png for example. 
            cover : The cover generated by the antenna_cover_greedy algorithm
            tol: The max allowable tolerance of uncovered points.
            coverage_units: The coverage units to convert the detection_coverage curves to (this is assumed to be given to us in km).
        """
        self.centers = [c[0] for c in cover   ]
        self.indexes = [c[1] for c in cover   ]
        self.state =  [c[2] for c in cover  ]
        # load the detection coverage file.
        self.detection_coverage = antennacover.read_detection_coverage(detection_coverage_file,coverage_units)
        # the polygons representing the antenna shapes (rotated and translated)
        self.cover_polygons = excessarea.generate_antenna_cover_polygons(self.indexes,self.state,self.centers,self.detection_coverage)
        # The bounding polygon consisting of the shore and the interference contour.
        self.bounding_polygon = bounding_polygon
        # Generate test points for computing area differences (polygon difference does not always work)
        self.generate_test_points()
        try:
            self.flag = False
            # BUGBUG -- this is to compensate for a bug in the data set.
            # Generate the test points to check for coverage.
            ratio = self.compute_uncovered_ratio(self.cover_polygons)
        except:
            # bounding polygon is not a simple polygon. Revert back to point counting.
            print "Exception -- fallback to slow routine"
            self.flag = True
            ratio =  self.compute_uncovered_ratio1(self.cover_polygons)
        # allow ourselves at least the given tolerance.
        self.max_ratio = max(ratio,tol)
        # The counter that allows us to step through the array of polygons repeatedly
        self.index_count = 0
        # The number of moves allowed.
        if steps == 0:
            self.steps = len(self.cover_polygons) * 50
        else:
            self.steps = steps

