# This software was developed by employees of the National Institute
# of Standards and Technology (NIST), an agency of the Federal
# Government. Pursuant to title 17 United States Code Section 105, works
# of NIST employees are not subject to copyright protection in the United
# States and are considered to be in the public domain. Permission to freely
# use, copy, modify, and distribute this software and its documentation
# without fee is hereby granted, provided that this notice and disclaimer
# of warranty appears in all copies.
#
# THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
# EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
# TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
# IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
# AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION
# WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE
# ERROR FREE. IN NO EVENT SHALL NASA BE LIABLE FOR ANY DAMAGES, INCLUDING,
# BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES,
# ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS
# SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR
# OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY
# OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT
# OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
#
# Distributions of NIST software should also include copyright and licensing
# statements of any third-party software that are legally bundled with
# the code in compliance with the conditions of those licenses.

import unittest
import sys
sys.path.append('../')
from circle import Circle
import line
from line import Line
import math
import numpy as np
import pdb

def areaOfTriangle(l1,l2,l3):
    a = l1.length()
    b = l2.length()
    c = l3.length()
    s = (a + b + c) /2
    area = math.sqrt(s*(s -a)*(s-b)*(s-c))
    return area

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

    def testLineStabsCircle2(self):
        p1 = [4,4]
        p2 = [10,8]
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 1)
        self.assertTrue(c is not None)


    def testLineIntersectsCircle(self):
        p1 = [0,6]
        p2 = [4,-1]
        from line import Line
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 2)
        self.assertTrue(c is not None)
        # Check if c endpoints are on the circle.
        # First translate to 0 coordinates
        c1 = np.subtract(c.get_p1(),self.circle.get_center())
        c2 = np.subtract(c.get_p2(),self.circle.get_center())
        r = math.sqrt(c1[0]**2 + c1[1]**2)
        self.assertTrue(np.isclose(r,self.circle.get_radius()))
        r = math.sqrt(c2[0]**2 + c2[1]**2)
        self.assertTrue(np.isclose(r,self.circle.get_radius()))

    def testLineIntersectsCircleComputeSlice(self):
        p1 = [0,6]
        p2 = [5,-1]
        from line import Line
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 2)
        self.assertTrue(c is not None)
        sliceArea = self.circle.compute_polar_slice_area([line_segment])
        self.assertTrue(sliceArea > 0)
        l1 = line.Line(self.circle.get_center(),c.get_p1())
        l2 = line.Line(self.circle.get_center(),c.get_p2())
        # The angle that is outside the wedge.
        angle =  2*math.pi -  l2.angle(l1)
        # The portion of the circle that is outside the wedge
        portion_outside = angle/(2*math.pi)*self.circle.area()
        # The triangle area
        triangleArea = areaOfTriangle(c,l1,l2)
        # The following is the exact computation of the excess area.
        segmentArea = self.circle.area() - (triangleArea + portion_outside)
        # Now check if the exact computation matches the approximate
        # computation approximately.
        self.assertTrue(np.allclose(segmentArea,sliceArea,rtol=1e-3))
        print("segmentArea " + str( segmentArea) +  " approximate sliceArea " 
                    + str(sliceArea) )


    def testLineIntersectsCircle2(self):
        p1 = [2,1]
        p2 = [9,7]
        line_segment = Line(p1,p2)
        b,l,c = self.circle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 2)
        self.assertTrue(c is not None)
        self.assertTrue(self.circle.on(c.get_p1()))
        self.assertTrue(self.circle.on(c.get_p2()))

    def testLineMissesCircle(self):
        p1 = [20,80]
        p2 = [40,70]
        line_segment = Line(p1,p2)
        self.circle = Circle(center=[60,40], radius=20)
        b,l,c = self.circle.collides(line_segment)
        self.assertFalse(b)
        self.assertTrue(len(l) == 1)
        self.assertTrue(c is None)

    #def testLineIntersectsCircle3(self):
    #    l2 = Line([50, 70], [58.78679656440357, 61.21320343559643])
    #    self.circle = Circle(center = [30, 60], radius = 22.360679775)
    #    b,l,c = self.circle.collides(l2)
    #    print "testLineIntersectsCircle3 : ", c

if __name__ == 'main':
    unittest.main()


