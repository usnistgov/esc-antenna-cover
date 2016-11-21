import unittest
import sys
sys.path.append('../')
import circle
import circles
import line

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

    
    def testMinimumCircleSetCoverForLineSet(self):
        line_endpoints = [[20,80],[40,70],[60,60],[60,40],[80,30],[90,20]]
        centers = [[10,70],[30,60],[50,55],[50,30],[60,20]]

        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))

        print "line_segments ", line_segments

        circles.min_area_cover_discrete(centers,line_segments)

    def testMinimumCircleSetCoverForLineSetGreedy(self):
        line_endpoints = [[20,80],[40,70],[60,60],[60,40],[80,30],[90,20]]
        centers = [[10,70],[30,60],[50,55],[50,30],[60,20]]

        line_segments = []
        for i in range( len(line_endpoints) - 1 ):
            line_segments.append(line.Line(line_endpoints[i],line_endpoints[i+1]))

        print "line_segments ", line_segments

        print circles.min_area_cover_greedy(centers,line_segments,[])

