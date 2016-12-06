import unittest
import sys
sys.path.append('../')
import circle
import circles
import line
import pdb
import random

def total_area(circle_collection):
    total_area = 0
    for c in circle_collection:
        total_area = total_area + c.get_radius()**2
    return total_area

def excess_area(circle_collection, covered_lines):
    total_area = 0
    for i in range(0,len(circle_collection)):
        total_area = total_area + circle_collection[i].compute_slice(covered_lines[i])
    return total_area

def printCover(lines,cover,centers,covered_segments,output_file):
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

    X = []
    Y = []
    r = []
    for c in centers:
        X.append(c[0])
        Y.append(c[1])
        r.append(1)

    f.write("X = " + str(X) + ";\n")
    f.write("Y = " + str(Y) + ";\n")
    f.write("r = 0.5*" + str(r) + ";\n")
    f.write("centers = [X',Y'];\n")
    f.write("viscircles (centers,r','Color','b');\n")

    # Draw the cover
    colors = ['b','r','g','y','k']
    for i in range(0,len(cover)):
        center = cover[i].get_center()
        radius = cover[i].get_radius()
        color = colors[i%len(colors)]
        f.write("viscircles (" + str(center) + "," + str(radius) + "," + "'Color'" + "," + "'" + str(color) + "');\n")
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
        r.append(1)
    f.write("r = 0.5*" + str(r) + ";\n")
    f.write("viscircles (centers,r','Color','g');\n")

    earea = excess_area(cover,covered_segments)
    f.write("total_area = " + str(total_area(cover)) + ";\n")
    print "total_area = " , str(total_area(cover))
    print ("excess_area = "+ str(earea))
    f.write( "excess_area = " +  str(earea) + ";\n")
    f.close()
    

class CirclesCoverTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCoversLine(self):
        l = line.Line([20,15],[37,37])
        c = [circle.Circle(center=[20,20], radius=10),circle.Circle(center=[35,25],radius=10),circle.Circle(center=[35,35],radius=5),circle.Circle(center=[20,50],radius=10)]
        b,r = circles.covers_line(c,l)
        self.assertTrue(b and len(r) == 3)

    def testCoversLine2(self):
        print "testCoversLine2"
        l = line.Line([15,45],[38,44])
        c = [ circle.Circle(center=[20,35], radius=13), 
                    circle.Circle(center=[35,40], radius=13) ]
        #c = [ circle.Circle(center=[20,35], radius=13), 
        #            circle.Circle(center=[35,40], radius=13), 
        #            circle.Circle(center=[45,55], radius=8),
        #            circle.Circle(center=[65,43], radius = 9),
        #            circle.Circle(center=[53,35], radius=10),
        #            circle.Circle(center=[40,25], radius=6) ]

        b,r = circles.covers_line(c,l)
        self.assertTrue(b)

    def testDoesNotCoverLine(self):
        l = line.Line([20,15],[37,37])
        c = [circle.Circle(center=[20,20], radius=10),circle.Circle(center=[35,25],radius=3),circle.Circle(center=[35,35],radius=5),circle.Circle(center=[20,50],radius=10)]
        b,r = circles.covers_line(c,l)
        self.assertFalse(b)

    def testDoesNotCoverLine2(self):
        l = line.Line([20,15],[57,57])
        c = [circle.Circle(center=[20,20], radius=10),circle.Circle(center=[35,25],radius=10),circle.Circle(center=[35,35],radius=5)]
        (b,r) = circles.covers_line(c,l)
        self.assertFalse(b)


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
        print "solnX = ",soln_x
        print "solnY = ",soln_y
        print "solnR = ",soln_r

    

    def testMinimumCircleSetCoverForLineSetGreedy(self):
        line_endpoints = [[20,80],[50,70],[60,60],[80,50],[60,40],[80,30],[100,30],[90,20]]
        centers = [[10,70],[30,60],[50,55],[50,30],[60,20]]
        savedCentrs = list(centers)
        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        circ,covered_segments = circles.min_area_cover_greedy(centers,line_segments)
        self.assertTrue(len(circ) == len(covered_segments))
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        for line_segment in line_segments:
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p1()):
                    covered = True
            self.assertTrue(covered)
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p2()):
                    covered = True
            self.assertTrue(covered)

        # Check if the segments in the cover are distinct subsets.
        for i in range(0,len(covered_segments)):
            for k in range(i+1,len(covered_segments) -1 ):
                for l in covered_segments[i]:
                    for m in covered_segments[k]:
                        self.assertFalse(l == m)
                

        printCover(line_segments,circ,savedCentrs,covered_segments,"testMinimumCircleSetCoverForLineSetGreedy.m")


    def testMinimumCircleSetCoverForLineSetGreedy2(self):
        line_endpoints = [[20,80],[50,70]]
        centers = [[10,70],[40,60]]
        savedCentrs = list(centers)
        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        circ,segments = circles.min_area_cover_greedy(centers,line_segments)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        for line_segment in line_segments:
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p1()):
                    covered = True
            self.assertTrue(covered)
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p2()):
                    covered = True
            self.assertTrue(covered)
        printCover(line_segments,circ,savedCentrs,segments,"testMinimumCircleSetCoverForLineSetGreedy2.m")

    def testMinimumCircleSetCoverForLineSetGreedy3(self):
        line_endpoints = [[20,80],[50,70]]
        centers = [[10,70],[30,60]]
        savedCentrs = list(centers)
        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        circ,segments = circles.min_area_cover_greedy(centers,line_segments)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        for line_segment in line_segments:
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p1()):
                    covered = True
            self.assertTrue(covered)
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p2()):
                    covered = True
            self.assertTrue(covered)
        printCover(line_segments,circ,savedCentrs,segments,"testMinimumCircleSetCoverForLineSetGreedy3.m")

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
            centers.append([p1,p2])
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        savedCentrs = list(centers)
        circ,included = circles.min_area_cover_greedy(centers,line_segments)
        for line_segment in line_segments:
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p1()):
                    covered = True
            self.assertTrue(covered)
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p2()):
                    covered = True
            self.assertTrue(covered)

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

        printCover(line_segments,circ,savedCentrs,included,"testMinimumCircleSetCoverForLineSetGreedyRandom.m")

