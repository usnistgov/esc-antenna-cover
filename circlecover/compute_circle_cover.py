
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

import circlecover
import argparse
import printcover


if __name__ == "__main__":
    # Read and parse the args.
    parser = argparse.ArgumentParser()
    parser.add_argument("-dist", type=int, default=0, help = "Min sensor spacing (km) default 0")
    parser.add_argument("-gs", type=int, default=100, help = "Grid size (default 100 changing this affects run time dramatically)")
    parser.add_argument("-pr", help="Definition of protected region units in meters",required=True)
    parser.add_argument("-of",default="output", help = "Output file name prefix")
    args = parser.parse_args()
    protection_region = args.pr
    min_ctr_dist = args.dist
    min_ctr_dist = args.dist
    output_file = args.of
    grid_size = args.gs

    # Load up the data.
    with open (protection_region, "r") as myfile:
        data=myfile.readlines()

    # Units are in Km. This should be converted to json format.
    esc_loc_x = [x/1000.0 for x in  eval(data[0])]
    esc_loc_y = [x/1000.0 for x in  eval(data[1])]
    ship_loc_x = [x/1000.0 for x in  eval(data[2])]
    ship_loc_y = [x/1000.0 for x in  eval(data[3])]

    possible_centers = []
    for i in range(0,len(esc_loc_x)):
        center = (esc_loc_x[i],esc_loc_y[i])
        possible_centers.append(center)

    interference_contour = []
    for i in range(0,len(ship_loc_x)):
        p = (ship_loc_x[i],ship_loc_y[i])
        interference_contour.append(p)

    testName = output_file
    #min_area_cover_greedy(possible_centers, interference_contour, min_center_distance=0):
    cover = circlecover.min_area_cover_greedy(possible_centers, interference_contour, min_center_distance=min_ctr_dist,ndivisions=grid_size)
    #printCover(interference_contour,cover,centers,min_separation,covered_segments,testName, algorithm):
    printcover.printCover(interference_contour, cover, possible_centers, min_ctr_dist,None,output_file,"AREA_COVER")

