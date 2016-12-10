import unittest
import sys
sys.path.append('../')
from line import Line
import numpy as np
import pdb

class LineTest(unittest.TestCase):

    def setUp(self):
        pass
    
    def testSplitLine(self):
        line = Line([5,4], [3,2])
        lines =  line.split(2)
        self.assertTrue(len(lines) == 2)
        s1 = lines[0].slope()
        s2 = lines[1].slope()
        self.assertTrue(np.isclose(s1,s2))

    def testIntersects(self):
        line1 = Line([5,4], [3,2])
        line2 = Line([7,8],[1,2])
        res, vals = line1.intersection(line2)
        self.assertFalse(res)
        self.assertTrue(vals is None)
        b = line1.intersects(line2)
        self.assertFalse(b)
	
    def testIntersects2(self):
        line1 = Line([5,4], [3,2])
        line2 = Line([3,4],[8,1])
        res, vals = line1.intersection(line2)
        self.assertTrue(res)
        self.assertTrue(vals is not None)
        self.assertTrue(vals == (4.25,3.25))
        self.assertTrue(line1.isCollinear(vals))

    def testIntersects4(self):
        line1 = Line([-1,-1], [1,1])
        line2 = Line([-1,1], [1,-1])
        res, vals = line1.intersection(line2)
        self.assertTrue(res)
        self.assertTrue(vals is not None)
        self.assertTrue(vals == (0,0))
        self.assertTrue(line1.isCollinear(vals))

    def testIntersects3(self):
        line1 = Line([5,4], [3,2])
        line2 = Line([3,4],[1,5])
        res, vals = line1.intersection(line2)
        self.assertFalse(res)
        self.assertTrue(vals is None)
        res, vals = line2.intersection(line1)
        self.assertFalse(res)
        self.assertTrue(vals is None)

    def testIsCollinear(self):
        line1 = Line([1,1],[3,3])
        self.assertTrue(line1.isCollinear([2,2]))

    def testIntersects5(self):
        #line1 = Line([13,6],[4,13])
        line1 = Line([4,13],[13,16])
        line2 = Line([4,3],[6,15])
        res1, vals1 = line1.intersection(line2)
        self.assertTrue(res1)
        self.assertTrue(vals1 is not None)
        res2, vals2 = line1.intersection(line2)
        self.assertTrue(res2)
        self.assertTrue(vals2 is not None)
        self.assertTrue(vals1 == vals2)
        self.assertTrue(line1.isCollinear(vals1))
        self.assertTrue(line2.isCollinear(vals1))


if __name__ == 'main':
    unittest.main()
