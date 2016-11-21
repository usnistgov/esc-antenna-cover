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
        c = [ circle.Circle(center=[25,35], radius=13), 
                    circle.Circle(center=[35,40], radius=13) ]
        #c = [ circle.Circle(center=[25,35], radius=13), 
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
        (b,r) = circles.min_area_cover_for_line(c,l)
        self.assertTrue(b)

    def testMinimumAreaCoverForLineSet(self):
        lines = [line.Line([15,45],[38,44]),line.Line([38,44],[58,38]),line.Line([58,38],[66,49])]
        cir = [ circle.Circle(center=[25,35], radius=13), 
                    circle.Circle(center=[35,40], radius=13), 
                    circle.Circle(center=[45,55], radius=8),
                    circle.Circle(center=[65,43], radius = 9),
                    circle.Circle(center=[53,35], radius = 8) ,
                    circle.Circle(center=[53,35], radius=10) ]
        b,c = circles.min_area_cover(cir,lines)
        self.assertTrue(b)
        

