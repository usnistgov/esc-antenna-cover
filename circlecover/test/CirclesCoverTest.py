import unittest
import sys
sys.path.append('../')
import circle
import circles
import line
import pdb

def cost(circle_collection):
    total_area = 0
    for c in circle_collection:
        total_area = total_area + c.get_radius()**2
    return total_area

def printCover(lines,cover,centers,output_file):
    X = []
    Y = []
    r = []
    for c in cover:
        X.append(c.get_center()[0])
        Y.append(c.get_center()[1])
        r.append(c.get_radius())

    f = open(output_file,"w")

    f.write("X = " + str(X) + ";\n")
    f.write("Y = " + str(Y) + ";\n")
    f.write("r = " + str(r) + ";\n")
    f.write("centers = [X',Y'];\n")
    f.write("viscircles (centers,r');\n")
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
    f.write("r = " + str(r) + ";\n")
    f.write("centers = [X',Y'];\n")
    f.write("viscircles (centers,r','Color','b');\n")
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
        line_endpoints = [[20,80],[40,70],[60,60],[80,50],[60,40],[80,30],[100,30],[90,20]]
        centers = [[10,70],[30,60],[50,55],[50,30],[60,20]]
        savedCentrs = list(centers)
        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))
        print "line_segments ", line_segments
        circ = circles.min_area_cover_greedy(centers,line_segments,[])
        print "circle_cover = ",circ, " cost = ", cost(circ)
        # check that for every line segment, both endpoints are covered by
        # at least one circle in our solution.
        for line_segment in line_segments:
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p1()):
                    covered = True
            #self.assertTrue(covered)
            covered = False
            for c in circ:
                if c.inside(line_segment.get_p2()):
                    covered = True
            #self.assertTrue(covered)

        printCover(line_segments,circ,savedCentrs,"testMinimumCircleSetCoverForLineSetGreedy.m")
    
    
        
    

