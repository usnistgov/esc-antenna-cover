
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
import circle
import circlecover
import line
import pdb
import random
import math
import logging
import numpy as np
from line import Line

def areaOfTriangle(l1,l2,l3):
    a = l1.length()
    b = l2.length()
    c = l3.length()
    s = (a + b + c) /2
    area = math.sqrt(s*(s -a)*(s-b)*(s-c))
    return area

class ExcessAreaTest(unittest.TestCase):

    def testCoverAreaNumericalIntegration(self):
        """
        Test the excess area computation.
        """
        p1 = [20,30]
        p2 = [50,60]
        line_segment = Line(p1,p2)
        lines = [line_segment]
        ccle = circle.Circle(center=[40,40], radius=15)
        b,l,c = ccle.collides(line_segment)
        self.assertTrue(b)
        self.assertTrue(len(l) == 2)
        self.assertTrue(c is not None)
        l1 = line.Line(ccle.get_center(),c.get_p1())
        l2 = line.Line(ccle.get_center(),c.get_p2())
     	angle = math.pi - l1.angle(l2)
        # The angle that is outside the wedge.
        remaining_angle =  2*math.pi -   angle
        # The portion of the circle that is outside the wedge
        portion_outside = remaining_angle/(2*math.pi)*ccle.area()
        # The triangle area
        triangleArea = areaOfTriangle(c,l1,l2)
        # The following is the exact computation of the excess area.
        segmentArea = ccle.area() - (triangleArea + portion_outside)
        # Now compare it with the numerical integration result.
        segmentArea1,totalarea = circlecover.compute_excess_area([ccle],lines)
        print "segmentArea ", segmentArea, " segmentArea1 ", segmentArea1
        self.assertTrue(np.allclose(segmentArea,segmentArea1,rtol=.01))
        self.assertTrue(np.allclose(ccle.area(),totalarea,rtol=.01))

    def testCoverAreaNumericalIntegration2(self):
        p1 = [99.0,29.0]
        p2 = [100.0,30.0]
        line_segment = line.Line(p1, p2)
        ccle = circle.Circle(center=[60,20], radius=41.231)
        self.assertTrue(ccle.inside(p1) and ccle.inside(p2))
        segmentArea,totalArea = circlecover.compute_excess_area([ccle],[line_segment],grid_divisions=200)
        print "segmentArea ", segmentArea
        self.assertTrue(segmentArea != 0)

