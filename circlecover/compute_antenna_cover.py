
import antennacover
import excessarea
import argparse
import simannealer
import printcover


if __name__ == "__main__":
    # Read and parse the args.
    parser = argparse.ArgumentParser()
    parser.add_argument("-dist", type=int, default=0, help = "Min sensor spacing (km) default 0")
    parser.add_argument("-gs", type=int, default=400, help = "Grid size (default 400)")
    parser.add_argument("-pr", help="Definition of protected region units in meters",required=True)
    parser.add_argument("-ap", help = "Definition of antenna patterns unit in Km.",required=True)
    parser.add_argument("-anneal", type = int, default=0, help="Number of steps to run the annealer")
    parser.add_argument("-of",default="output", help = "Output file name prefix")
    parser.add_argument("-to",type=float,default=.005, help = "outage tolerance (default = .005)")
    args = parser.parse_args()
    protection_region = args.pr
    do_anneal = args.anneal != 0
    coverage_file = args.ap
    #min_ctr_dist = args.dist
    min_ctr_dist = args.dist
    output_file = args.of
    grid_size = args.gs
    tol = args.to
    antennacover.NDIVISIONS=grid_size
    
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

    bounding_polygon = excess_area.generate_bounding_polygon(possible_centers,interference_contour)
    testName = output_file
    cover = antennacover.min_antenna_area_cover_greedy(possible_centers, bounding_polygon, coverage_file, min_center_distance=min_ctr_dist,tol=tol)
    printcover.printAntennaCover(output_file,bounding_polygon, possible_centers, cover,coverage_file,min_ctr_dist)
    if do_anneal:
        annealr = simannealer.SimAnneal(bounding_polygon, possible_centers, coverage_file,cover,steps=args.anneal,tol=tol)
        annealr.anneal()
        cover = annealr.get_result()
        printcover.printAntennaCover(output_file + "Anneal", bounding_polygon, possible_centers, cover,coverage_file,min_ctr_dist)


