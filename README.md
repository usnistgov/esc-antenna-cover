## Geometric Circle-cover

Minimum area circle cover from a set of center points covering a set of lines. 

Given a set of M points on a plane where circles can be centered and
a set of N line segments which need to be covered by the circles,
find the minimum area circle cover for the line segments. That is,
find the radii of the circles and the centers (chosen from the M points)
such that all the N line segments are covered and the total area of the
circles is minimized.

Note that a line segment is covered if no part of it is OUTSIDE a circle.

The python code in this project implements the following algorithm:

     1. Find the worst line, i.e. the line requiring the largest additional
     circle area for its best circle center option with the corresponding
     line segment entirely in the circle. 

     2. Construct / extend the corresponding circle. 

     3. Remove fully covered line segments from the set and cut partially 
     covered segments at the circle boundary. Remove the center from the set of possible centers. 
     Iterate until no more line segments remain.



## Installing the code

Install Shapley Package 

For Linux:

    pip install shapely


For windows: Download the shapley installer for your architecture:

    http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely

Run pip install on the downloaded .whl file.

Now run setup for the current package:

    python setup.py install

## USAGE
### Python

### MATLAB

Add the circlecover directory to your matlab path.

See the documentation help min\_cover for usage.

See the example in circlecover/test/CircleCoverTest.m


    esc_loc_x = [1771380,1769310,1769790,1768380,176739,1764690,...
        1762020,1759920,1753110,1741950,1752210,1757010,1761870,...
        1768230,1772820,1777110,1781610,1786920,1793220];

    esc_loc_y = [1827030,1817070,1806990,1797090,1787100,1776840,...
        1767270,1756950,1746690,1735050,1727220,1717290,1707360,...
        1697370,1687320,1677450,1667400,1657350,1647360];

    ic_x = [1847012,1844913,1845660,1834150,1823280,1811715,...
        1807512,1806671,1810710,1807769,1817910,1822503,1827218,...
        1823623,1828432,1842183,1846928,1852378,1858591];

    ic_y = [1843636,1833617,1823583,1811442,1799284,1787072,1777140,...
        1767066,1759078,1749183,1741311,1731358,1721401,1709309,...
        1699318,1691518,1681523,1671542,1661589];


    esc_loc = [esc_loc_x',esc_loc_y'];

    ic = [ic_x',ic_y'];

    distance = 60;

    [centers_x,centers_y,radius] = min_cover(esc_loc,ic,distance);

    disp('centers_x');
    disp(centers_x);
    disp('centers_y');
    disp(centers_y);
    disp('radius');
    disp(radius);




## Disclaimers

This software was developed by employees of the National Institute
of Standards and Technology (NIST), an agency of the Federal
Government. Pursuant to title 17 United States Code Section 105, works
of NIST employees are not subject to copyright protection in the United
States and are considered to be in the public domain. Permission to freely
use, copy, modify, and distribute this software and its documentation
without fee is hereby granted, provided that this notice and disclaimer
of warranty appears in all copies.

THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION
WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE
ERROR FREE. IN NO EVENT SHALL NASA BE LIABLE FOR ANY DAMAGES, INCLUDING,
BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES,
ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS
SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR
OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY
OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT
OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.

Distributions of NIST software should also include copyright and licensing
statements of any third-party software that are legally bundled with
the code in compliance with the conditions of those licenses.



Acknowledgement
===============


-- The algorithm implemented here was suggested by Stefan Haustein see:

	http://stackoverflow.com/questions/40748412/minimun-area-geometric-cover-for-a-set-of-line-segments


-- The problem of minimum excess area circle cover was proposed by Tim Hall. This was modified to minimum area cover
   as an approximation.

-- Improvements to the algorithm were evolved and bad ideas pruned as a result of discussions with Anastase Nakassis at NIST.

-- This algorithm was developed for placement of and clibration of ESC sensors for 3.5 GHz spectrum sharing
   but it is generally applicable for geometric line cover.

-- The Author thanks Thao Nguyen, Tim Hall and Anirudha Sahoo for proposing the ESC Sensor placement problem.

Future Work
==========

Worst case sub-optimality analysis.

