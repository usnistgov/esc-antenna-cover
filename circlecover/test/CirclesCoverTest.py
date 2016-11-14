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
        l = line.line([20,15],[37,37])
        c = [circle.circle(center=[20,20], radius=10),circle.circle(center=[35,25],radius=10),circle.circle(center=[35,35],radius=5),circle.circle(center=[20,50],radius=10)]
        self.assertTrue(circles.covers_line(c,l))

    def testDoesNotCoverLine(self):
        l = line.line([20,15],[37,37])
        c = [circle.circle(center=[20,20], radius=10),circle.circle(center=[35,25],radius=3),circle.circle(center=[35,35],radius=5),circle.circle(center=[20,50],radius=10)]
        self.assertFalse(circles.covers_line(c,l))

    def testDoesNotCoverLine2(self):
        l = line.line([20,15],[57,57])
        c = [circle.circle(center=[20,20], radius=10),circle.circle(center=[35,25],radius=10),circle.circle(center=[35,35],radius=5)]
        self.assertFalse(circles.covers_line(c,l))

