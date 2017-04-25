
import circlecover
import argparse
import printcover


if __name__ == "__main__":
    # Read and parse the args.
    parser = argparse.ArgumentParser()
    parser.add_argument("-dist", type=int, default=0, help = "Min sensor spacing (km) default 0")
    parser.add_argument("-gs", type=int, default=400, help = "Grid size (default 400)")
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
    cover,covered = circlecover.min_area_cover_greedy(possible_centers, interference_contour, min_center_distance=min_ctr_dist,ndivisions=grid_size)
    #printCover(interference_contour,cover,centers,min_separation,covered_segments,testName, algorithm):
    printcover.printCover(interference_contour, cover, possible_centers, min_ctr_dist,None,output_file,"AREA_COVER")

