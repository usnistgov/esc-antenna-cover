import unittest
import sys
sys.path.append('../')
from line import Line
import numpy as np

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
	
    def testIntersects2(self):
        line1 = Line([5,4], [3,2])
        line2 = Line([3,4],[8,1])
        res, vals = line1.intersection(line2)
        self.assertTrue(res)
        self.assertTrue(vals is not None)
        self.assertTrue(vals == [4,3])

    def testIntersects4(self):
        line1 = Line([-1,-1], [1,1])
        line2 = Line([-1,1], [1,-1])
        res, vals = line1.intersection(line2)
        self.assertTrue(res)
        self.assertTrue(vals is not None)
        self.assertTrue(vals == [0,0])

    def testIntersects3(self):
        line1 = Line([5,4], [3,2])
        line2 = Line([3,4],[1,5])
        res, vals = line1.intersection(line2)
        self.assertFalse(res)
        self.assertTrue(vals is None)
        res, vals = line2.intersection(line1)
        self.assertFalse(res)
        self.assertTrue(vals is None)

if __name__ == 'main':
    unittest.main()
