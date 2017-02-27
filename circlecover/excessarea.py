from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
import circle
import line
import antennacover
import copy
import pdb
from collections import namedtuple

def generate_bounding_polygon(possible_centers,interference_contour):
    """ generate a bounding polygon consisting of possible_centers and interference_contour points. """
    centers = copy.copy(possible_centers)
    # Create a multipoint polygon coonsisting of the original contour 
    # and the possible centers
    points = [point for point in interference_contour]
    # The centers and the shore points are listed in the same sorted order
    centers.reverse()
    for point in centers:
        points.append(point)
    # The polygon encloses the interference contour as well as the shore.
    # This is the region that needs to be covered.
    mp = Polygon(points)
    return mp

def compute_excess_area_for_antenna_cover(indexes, angles, centers, detection_coverage_file,
            possible_centers, interference_contour):
    """
    Parameters:
        indexes - the selected contours from detection_coverage_file
        angles - the orientation of each antenna.
        centers : the centers of the sensors (ordered).
        interference_contour : The interference contour point set (ordered)
        detection_coverage_file: The detection coverage file.
    """



    def generate_antenna_cover_polygons(indexes,angles,centres,detection_coverage):
        polygons = []
        for i in range(0,len(indexes)):
            polygon = detection_coverage[indexes[i]].lobe
            angle = angles[i]
            center = centers[i]
            rotated_translated_polygon = antennacover.translate(antennacover.rotate(polygon,angle),center)
            polygons.append(rotated_translated_polygon)
        return polygons

    def point_covered_by_lobe(point,antenna_cover_lobes):
        for i in range(0,len(antenna_cover_lobes)):
            if antenna_cover_lobes[i].contains(point):
                return True
        return False
            




    # load the detection coverage file.
    detection_coverage = antennacover.read_detection_coverage(detection_coverage_file)
    
    # the bounding polygon representing the if-contour and the possible centers. This bounds the region that
    # needs to be covered.
    bounding_polygon = generate_bounding_polygon(possible_centers,interference_contour)

    # Line string representing the interference contour. Add first and last point 
    # of possible_centers to this. This is the interference contour that is on the sea.
    interference_contour.append(possible_centers[-1])
    interference_contour.insert(0,possible_centers[0])
    interference_contour_linestring = LineString(interference_contour)
    
    # Line string representing the possible sensor locations 
    possible_centers_linestring = LineString(possible_centers)

    # the polygons representing the antenna shapes (rotated and translated)
    antenna_cover_polygons = generate_antenna_cover_polygons(indexes,angles,centers,detection_coverage)
    # Take the union of these shapes
    union = antenna_cover_polygons[0]
    for i in range(1,len(antenna_cover_polygons)):
        union = union.union(antenna_cover_polygons[i])
    minx,miny,maxx,maxy = union.bounds

    # Generate a point set and classify.
    ndivs = 100
    deltax,deltay = float(maxx-minx)/ndivs, float(maxy-miny)/ndivs
    area_per_grid_point = deltax*deltay
    
    excess_sea_coverage_count = 0
    excess_land_coverage_count = 0
    for i in range(0,ndivs):
        for j in range(0,ndivs):
            p = Point(minx + i*deltax, miny + j*deltay)
            if point_covered_by_lobe(p,antenna_cover_polygons) and not bounding_polygon.contains(p):
                # The point p is now either on sea or land.
                d1 = interference_contour_linestring.distance(p)
                d2 = possible_centers_linestring.distance(p)
                if d1 <  d2 :
                    excess_sea_coverage_count = excess_sea_coverage_count + 1
                else:
                    excess_land_coverage_count = excess_land_coverage_count + 1

    return (excess_sea_coverage_count*area_per_grid_point, excess_land_coverage_count*area_per_grid_point)
                        
    

        
        
        
        



def compute_excess_area(circles, line_segments, grid_divisions=200):
    """
    This is a support function that is used to evaluate a cover algorithm.
    compute the excess area - i.e. the area between the line collection
    and the perephry of the circle using numerical integration.
    Returns a tuple -- the area and the grid size used for computation.
    The grid size is the length of each grid square.
    """

    def isoutside(circle, segments,point):
        """
        return true if the point is outside the line segment (between the
        lines and the circumference of the circle.
        """
        l = line.Line(circle.get_center(),point)
        intersection_count = 0
        for line_segment in segments:
            if  line_segment.intersects(l):
                intersection_count = intersection_count + 1
        #odd number of intersections means point is in 
        #excess area region.
        # Note that this ONLY works if the interference contour
        # is entirely to ONE side of the centers (as will usually 
        # be the case
        return (intersection_count % 2) == 1

    def get_line_segments_for_circle(circle,line_segments):
        """
        Compute the intersection of the circle with a set of line segments.
        """
        retval = []
        for l in line_segments:
            b,excluded,included = circle.collides(l)
            if b:
                retval.append(included)
        return retval


    def generate_grid(circles, grid_divisions):
        """
        Generate a grid for numerical integration.
        Returns a list of points enclosed in the circles.
        """
        centers = []
        left_bottom = None
        right_top = None
        for c in circles:
            if left_bottom == None:
                left_bottom = [float(c.get_center()[0] - c.get_radius()), float(c.get_center()[1] - c.get_radius())]
            else:
                if c.get_center()[0] - c.get_radius() < left_bottom[0]:
                    left_bottom[0] = c.get_center()[0] - c.get_radius()
                if c.get_center()[1] - c.get_radius() < left_bottom[1]:
                    left_bottom[1] = c.get_center()[1] - c.get_radius()

            if right_top == None:
                right_top = [float(c.get_center()[0] + c.get_radius()), float(c.get_center()[1] + c.get_radius())]
            else:
                if c.get_center()[0] + c.get_radius() > right_top[0]:
                    right_top[0] = c.get_center()[0] + c.get_radius()
                if c.get_center()[1] + c.get_radius() > right_top[1]:
                    right_top[1] = c.get_center()[1] + c.get_radius()

        # find the side of the grid. We insist that each dimension should be broken up into
        # atleast grid_intervals intervals (arbitrary).
        x_distance = right_top[0] - left_bottom[0]
        y_distance = right_top[1] - left_bottom[1]

        if y_distance > x_distance:
            grid_size = y_distance/float(grid_divisions)
        else:
            grid_size = x_distance/float(grid_divisions)


        x_divisions = int(x_distance/grid_size)
        y_divisions = int(y_distance/grid_size)

        grid = []
        for i in range (0,x_divisions):
            for j in range (0,y_divisions):
                point = [left_bottom[0] + i*grid_size,left_bottom[1] + j*grid_size]
                grid.append(point)

        # Find the subset of the grid that is inside at least one circle.
        newgrid = []
        for point in grid:
            found = False
            for c in circles:
                #if we find a circle that includes the point, 
                # then incude in our integration side.
                if not found and c.inside(point) :
                    newgrid.append(point)
                    found = True

        # newgrid is a list of points within the circle set.
        # we can integrate over this grid. grid_size is the dimension of one side
        # of the integration square.
        return newgrid,grid_size

    grid,grid_size = generate_grid(circles,grid_divisions)
    grid_area = len(grid)*grid_size*grid_size

    count = 0
    for c in circles:
        segments = get_line_segments_for_circle(c,line_segments)
        # figure out the points in the annulus. Note that we consider a point to be in the annulus if 
        # the center is in corresponding grid point is in the annulus
        filteredPoints = filter(lambda p: c.inside((p[0] + grid_size/2, p[1] + grid_size/2)) and isoutside(c,segments,(p[0] + grid_size/2, p[1] + grid_size/2)),grid)
        k = 0
        for p in filteredPoints:
            grid.remove(p)
            k = k + 1
        count = count + k
    area = count * grid_size * grid_size
    return area,grid_area
