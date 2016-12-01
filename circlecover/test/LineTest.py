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
	print lines
	print s1,s2
	self.assertTrue(np.isclose(s1,s2))
	

if __name__ == 'main':
    unittest.main()
