import circle
import line
import itertools
import sys
import pdb
import numpy as np
import traceback
import math

def covers_line(circles,l):
    """
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
        t,lines =  top.collides(l)
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
    #pdb.set_trace()
    return covers_line_worker(circles,l,[])


def cost(circle_list):
    """
    return the cost of a circle (its area)
    """
    c = 0
    for c in circle_list:
        c = c.r**2
    return c

def min_area_cover_for_line_brute_force(circles,l):
    """
    The minimum area cover of a line selected from a set of circles.
    This is truly a brutish brute force algorithm
    """
    permuted_circles = itertools.permutations(circles)
    currentMinCover = None
    currentMinCost =  sys.maxint
    found = False
    for t in permuted_circles:
        found,cover =  covers_line(t,l)
        #if found:
        #   print "cover ", cover, " cost " , cost(cover), " min_cost ", currentMinCost, " current Min Cover ", currentMinCover
        if found and cost(cover) < currentMinCost:
           currentMinCover = cover
           currentMinCost = cost(cover)
    return found,currentMinCover

def min_area_cover_for_line_greedy(circles,l):
    """
    The minimum area cover of a line selected from a set of circles.
    This is truly a brutish brute force algorithm
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
                                cover.remove(c)
                            coverset.remove(r)
                            coverset.append(c)
            coverset = coverset + cover
    return True, list(set(coverset))
        
    
def min_area_cover_greedy(centers,lines,cover):
    """

    See http://stackoverflow.com/questions/40748412/minimun-area-geometric-cover-for-a-set-of-line-segments

    Find the worst point, i.e. the point requiring the largest additional circle area for its best circle center 
    option with the corresponding line segment at least partially in the circle. 
    Construct / extend the corresponding circle.
    Remove fully covered line segments from the set and cut partially covered segments at the circle boundary.
    Iterate until no more line segments remain.
    """

    def distance(p1,p2):
        return math.sqrt((p1[0]  - p2[0])**2 + (p1[1] - p2[1])**2)
        
    def intersects_lines(c,line_set):
        retval = []
        for line in line_set:
            b,l = c.collides(line)
            retval = retval + l
        return retval

    def cost(circle_collection):
        total_area = 0
        for c in circle_collection:
            total_area = total_area + c.get_radius()**2
        return total_area

    max_min_distance =  -1
    max_min_center = None


    for line in lines:
        min_distance = 1e6
        mincenter = None
        # find the cheapest cover for each line. 
        for center in centers:
            dist = max(distance(center,line.get_p2()), distance(center,line.get_p1()))
            if dist < min_distance:
                min_distance = dist
                mincenter = center
        if min_distance > max_min_distance:
           max_min_distance = min_distance
           max_min_center = mincenter


    # At this point we have computed the max_min_distance over all lines and centers.
    c = circle.Circle(max_min_center,max_min_distance)
    cover.append(c)
    newlines = intersects_lines(c,lines)
    if len(newlines) == 0 :
        print "cost = ",cost(cover)
        return cover
    else:
        if len(centers) != 1:
            centers.remove(max_min_center)
        return min_area_cover_greedy(centers,newlines,cover)
    
                


def min_area_cover_discrete(centers,lines,brute=False):
    """
    The minimum area cover for a set of lines 
    from circles centered at given locations.
    This is a brute force algorithm that tries to find
    the optimum radius assignments.
    
    Parameters:
    
        - centers: A list of coordinates for the centers of the sensor
        - lines : A list of line segments denoting the protection zone.
        - brute: A flag that tells us whether to use brute force search or not.

    Returns:
        
        - A list of circles with minimum total area that covers the lines.

    """
    def distance(p1,p2):
        return math.sqrt((p1[0]  - p2[0])**2 + (p1[1] - p2[1])**2)
        
    def cost(circle_collection):
        total_area = 0
        for c in circle_collection:
            total_area = total_area + c.get_radius()**2
        return total_area

    max_distance = -1
    min_distance = 1e6
    # find the max distance from any center 
    for line in lines:
        for center in centers:
            d = distance(center,line.get_p1())
            if d > max_distance:
                max_distance = d
            if d < min_distance:
                min_distance = d
            d = distance(center,line.get_p2())
            if d > max_distance:
                max_distance = d
            if d < min_distance:
                min_distance = d

    delta_r = (max_distance - min_distance) / 5
    
    # divide the range up into 5 segments (corresponding to 5 sensor powers)
    radii  = np.arange(min_distance,max_distance + delta_r, delta_r)
    print "radii ", radii

    # Generate a permutation matrix of all possible radius permutations.
    radius_assignments = [p for p in itertools.product(radii,repeat=len(centers))]
    
    min_cost = 1e6
    min_cover = None

    try:
        for radius_assignment in radius_assignments:
            for i in range(0,len(centers)):
                circles_r = []
                for r in radius_assignment:
                    circle_i = circle.Circle(center = centers[i], radius=r)
                    circles_r.append(circle_i)
                if brute:
                    t,cover = min_area_cover_for_lines_brute_force(circles_r,lines)
                else:
                    t,cover = min_area_cover_for_lines_greedy(circles_r,lines)
                if t:
                    try:
                        _cost = cost(cover)
                    except:
                        pdb.set_trace()
                    if _cost < min_cost:
                        min_cost =  _cost
                        min_cover = cover

        print "min_cover " , str(min_cover), " min_cost ", str(min_cost)

        return min_cover
    except:
        traceback.print_exc()
        raise
        
    
            

