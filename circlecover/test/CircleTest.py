import unittest
import sys
sys.path.append('../')
from circle import Circle
from line import Line

class CircleTest(unittest.TestCase):

    def setUp(self):
        self.circle = Circle(center=[5,4], radius=4)

    def testLineInsideCircle(self):
        p1 = [3,4]
        p2 = [4,4]
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 0)
        self.assertTrue(c is not None)

    def testLineStabsCircle(self):
        p1 = [0,6]
        p2 = [4,4]
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 1)
        self.assertTrue(c is not None)
        print '-------------------------------------'
        for li in l :
            print str(li)
        print c

    def testLineStabsCircle2(self):
        p1 = [4,4]
        p2 = [10,8]
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 1)
        self.assertTrue(c is not None)
        print '-------------------------------------'
        for li in l :
            print str(li)
        print c


    def testLineIntersectsCircle(self):
        p1 = [0,6]
        p2 = [4,-1]
        from line import Line
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 2)
        self.assertTrue(c is not None)
        print '-------------------------------------'
        for li in l :
            print str(li)
        print c

    def testLineIntersectsCircle2(self):
        p1 = [2,1]
        p2 = [9,7]
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 2)
        self.assertTrue(c is not None)
        print '-------------------------------------'
        for li in l :
            print str(li)
        print c

    def testLineMissesCircle(self):
        p1 = [20,80]
        p2 = [40,70]
        line_segment = Line(p1,p2)
        self.circle = Circle(center=[60,40], radius=20)
        b,l,c = self.circle.collides(line_segment)
        self.assertFalse(b)
        self.assertTrue(len(l) == 1)
        self.assertTrue(c is None)

if __name__ == 'main':
    unittest.main()


