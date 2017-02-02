
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
