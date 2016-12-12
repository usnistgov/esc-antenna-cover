import circle as ccle
import line
import itertools
import sys
import pdb
import numpy as np
import traceback
import math
import random
from line import Line
import logging
import copy

logger = logging.getLogger("circlecover")



def covers_line(circles,l):
    """
    We are given a set of circles and ask for the subset of the
    circles that covers the line l.
    
    return true if the set of circles covers the line.
  
    Parameters:
        - circles the set of circles to test.
        - cover - the subset of the circles that covers l
    
    Returns:

        - a pair of <true,array[circle]> representing the cover that the
        algorithm found or  <false,None> if no such cover exists.

    """

    def covers_line_worker(circles,l,cover):
        if len(circles) ==  0:
            return False, None
        # Get the topmost circle on the list
        top = circles[0]
        # check if it collides with the given line sebment.
        # lines contains the pieces caused by the collision.
        t,lines,contained =  top.collides(l)
        # Construct a new list of circles by removing the topmost element of the list.
        newc = list(circles)
        newc.remove(top)
       
        #the circle did not collide with the segment.
        if not t:
            # Test if the rest of the circles cover the line segment
            return covers_line_worker(newc,l,cover)

        # If the circle did collide with the line segment, include it in the cover.
        cover.append(top)
        # Now check if there are remaining line segments. If there are none, the line
        # is completely enclosed .
        if len(lines) == 0:
            # completely enclosed.
            return True,cover
        if len(lines) == 1:
            # If there is one segment, get the cover for that segment
            remainder = lines[0]
            return covers_line_worker(newc,remainder,cover)
        elif len(lines) == 2:
            # If there are two one segments, get the cover for those segments.
            left = lines[0]
            right = lines[1]
            b1, newCover =  covers_line_worker(newc,left,cover) 
            if not b1:
                return False,None
            else:
                return covers_line_worker(newc,right,newCover)
    return covers_line_worker(circles,l,[])



def min_area_cover_for_line_brute_force(circles,l):
    """
    The minimum area cover of a line selected from a set of circles.
    This is a brutish brute force algorithm. It just tries to cover
    the lines using a set of circles by permuting the order of covering
    using different orders.

    Here we are given a set of circles and want to find the minimum cost cover
    for a line using brute force search.
    """
    def cost(circle_list):
        """
        return the cost of a circle (its area)
        """
        cst = 0
        for c in circle_list:
            cst = cst +  c.r**2
        return cst

    permuted_circles = itertools.permutations(circles)
    currentMinCover = None
    currentMinCost =  sys.maxint
    found = False
    for t in permuted_circles:
        found,cover =  covers_line(t,l)
        if found and cost(cover) < currentMinCost:
           currentMinCover = cover
           currentMinCost = cost(cover)
    return found,currentMinCover

def min_area_cover_for_line_greedy(circles,l):
    """
    The minimum area cover of a line selected from a set of circles.

    Sort the circles by radius and cover with the smallest radius first.
    """
    circles.sort(key=lambda x : x.get_radius)
    return covers_line(circles,l)

def min_area_cover_for_lines_brute_force(circles,lines):
    """
    The minimum area cover for a set of lines given
    a set of circles.
    """
    # coverset is the set of circles that covers the
    # given line set.
    coverset = []
    for line in lines:
        found,cover = min_area_cover_for_line_brute_force(circles,line)
        if not found:
            return False,None
        else:
            for c in cover:
                for r in coverset:
                    if r.get_center() == c.get_center():
                        if r.get_radius() < c.get_radius():
                            # Pick the larger circle and add it to the
                            # cover set.
                            if c in cover:
                                cover.remove(c)
                            coverset.remove(r)
                            coverset.append(c)
            coverset = coverset + cover
    return True, list(set(coverset))

        
def min_area_cover_for_lines_greedy(circles,lines):
    """
    The minimum area cover for a set of lines given
    a set of circles.
    The approach is to just cover each line in the lines set
    using the given circles starting from the smallest one.
    Then, take the uninon of these covers.

    This is for the situation where you know the circles and their
    radii apriori and you want to return the minimum area cover
    The algorithm proceeds by starting with the smallest circle and
    covering the lines, then recursively apply to the remaining lines.
    
    """
    # coverset is the set of circles that covers the
    # given line set.
    coverset = []
    for line in lines:
        found,cover = min_area_cover_for_line_greedy(circles,line)
        if not found:
            return False,None
        else:
            for c in cover:
                for r in coverset:
                    if r.get_center() == c.get_center():
                        if r.get_radius() < c.get_radius():
                            # Pick the larger circle and add it to the
                            # cover set.
                            if c in cover:
                                # Need to check here to ensure
                                # c is stil in cover - otherwise remove fails.
                                cover.remove(c)
                            coverset.remove(r)
                            coverset.append(c)
            coverset = coverset + cover
    return True, list(set(coverset))
        

def compute_excess_area(circles, line_segments):
    """
    compute the excess area - i.e. the area between the line collection
    and the perephry of the circle using numerical integration.
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


    def generate_grid(circles):
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
        # atleast 200 intervals (arbitrary).
        x_distance = right_top[0] - left_bottom[0]
        y_distance = right_top[1] - left_bottom[1]

        if y_distance > x_distance:
            grid_size = y_distance/200.0
        else:
            grid_size = x_distance/200.0


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

    grid,grid_size = generate_grid(circles)
    grid_area = len(grid)*grid_size*grid_size

    count = 0
    for c in circles:
        segments = get_line_segments_for_circle(c,line_segments)
        filteredPoints = filter(lambda p: c.inside(p) and isoutside(c,segments,p),grid)
        for p in filteredPoints:
            grid.remove(p)
        count = count + len(filteredPoints)
    logger.debug("grid_area : " + str(len(grid)*grid_size*grid_size))
    area = count * grid_size * grid_size
    return area,grid_area
        

        


def min_area_cover_greedy(possible_centers,line_segments, min_center_distance = 0, nsegments=10):
    """
    Greedy minimum area cover.

    See discussion on 

    http://stackoverflow.com/questions/40748412/minimun-area-geometric-cover-for-a-set-of-line-segments

    Find the worst line, i.e. the line requiring the largest additional
    circle area for its best circle center option with the corresponding
    line segment entirely in the circle.

    Construct / extend the corresponding circle.

    Remove fully covered line segments from the set and cut partially
    covered segments at the circle boundary.  

    Remove the center from the set of possible centers.
    
    Iterate until no more line segments remain.


    Parameters:

    possible_centers : a list of points where the circles can be placed.
    line_segments: a list of Line which need to be optimally covered.
    min_center_distance: distance constraint for minimum distance between centers (defaults to 0).
    nsegments : # of segments to break up each line (defaullts to 10).
    

    Returns:

    return a tuple consisting of two arrays

    -  list of circles that completely covers the given lines.
    -  a collection of line segments for the lines enclosed by the corresponding 
       circles.

    # TODO - remove the nsegments parameter
    

    """
    def distance(p1,p2):
        return math.sqrt((p1[0]  - p2[0])**2 + (p1[1] - p2[1])**2)
        
    def intersects_lines(c,line_set):
        """
        compute intersection of line set with circle and return a lists
        of segments that are excluded and included in the circles.
        """
        retval = []
        included_set = []
        for line in line_set:
            b,l,included = c.collides(line)
            retval = retval + l
            if b:
                included_set.append(included)
        return retval,included_set


    def improve_solution(cover,lines):
        """
        Given set of circes and a segment cover, check if cover can be improved
        by moving segments around between circles. We check the intersection areas
        between circles and ask if we can improve the solution by moving the segments
        in the intersection area to the smaller of the two circles.

        """
        print "Improve solution "
        circles = [ccle.Circle(center=c.get_center(), radius=c.get_radius()) for c in cover]
        covered = [ [line.Line(l.get_p1(),l.get_p2()) for l in li] for li in lines]
        # Check if we can find a better fit.
        # Move lines into smaller circles if possible.
        for i in range(0,len(circles)-1):
            for j in range(i+1,len(circles)):
                # a set of indices for lines that are enclosed in a smaller dia circle.
                try:
                    inds = [k for (k,c) in enumerate(covered[i]) if circles[i].overlaps(circles[j]) and  circles[j].encloses(c)]
                    if len(inds) != 0:
                        # Enhance the covered set for the j'th circle.
                        covered[j] = covered[j] + [covered[i][k] for k in inds]
                        # resize the j'th circle.
                        center = circles[j].get_center()
                        # Enhance the covered set for the j'th circle.
                        covered[j] = covered[j] + [covered[i][k] for k in inds]
                        # resize the j'th circle.
                        center = circles[j].get_center()
                        radius = np.max([max(distance(center,l.get_p1()),distance(center,l.get_p2())) for l in covered[j]])
                        circles[j].set_radius(radius)
                        for index in sorted(inds, reverse=True):
                            del covered[i][index]
                        center = circles[i].get_center()
                        if len(covered[i]) == 0:
                            radius = 0
                        else:
                            radius = np.max([max(distance(center,l.get_p1()),distance(center,l.get_p2())) for l in covered[i]])
                        circles[i].set_radius(radius)
                except:
                    pdb.set_trace()
        inds = [k for (k,c) in enumerate(circles) if c.get_radius() == 0]
        for i in sorted(inds, reverse=True):
            del covered[i]
            del circles[i]
        return circles,covered
    
    def min_area_cover_greedy_worker(centers,lines,cover,segments,disance_constraint):

        max_min_distance =  -1e6
        max_min_center = None


        for line in lines:
            min_distance = 1e6
            mincenter = None
            # find the cheapest cover for each line. 
            for center in centers:
                # The circle that can cover both endpoints of the line.
                dist = max(distance(center,line.get_p2()), distance(center,line.get_p1()))
                if dist < min_distance:
                    min_distance = dist
                    mincenter = center
            # max_min_distance is the biggest minimum circle that covers any line in the
            # collection.
            if min_distance > max_min_distance:
                max_min_distance = min_distance
                max_min_center = mincenter


        c = ccle.Circle(max_min_center,max_min_distance)

        # now compute how much of lines there is left behind and how much is
        # included in the cover.

        newlines,included_lines = intersects_lines(c,lines)

        cover.append(c)
        segments.append(included_lines)
    


        if len(newlines) == 0 :
            # no more line segments left to cover -- we are done.
            return cover,segments
        else:
            # If there are remaining line segments, iterate on the portion 
            # that is left.
            centers.remove(max_min_center)
            # Remove all the centers within our distance constraint so our constraint
            # is satisfied.
            centers_to_remove = [k for k,cntr in enumerate(centers) if distance(cntr,c.get_center()) < min_center_distance ]
            for k in sorted(centers_to_remove, reverse=True):
                del centers[k]
            return min_area_cover_greedy_worker(centers,newlines,cover,segments,min_center_distance)

    cover = [] 
    included_segments = []
    split_line_segments = []
    for l in line_segments:
         #segments = l.split(random.randint(1,nsegments))
         segments = l.split(nsegments)
         for s in segments:
             split_line_segments.append(s)


    cover,covered = min_area_cover_greedy_worker(possible_centers,split_line_segments,cover,included_segments,min_center_distance) 

    cover,covered = improve_solution(cover,covered)
    
    return cover,covered
