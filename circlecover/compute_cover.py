
import antennacover
import excessarea
import argparse
import simannealer
import printcover


if __name__ == "__main__":
    # Read and parse the args.
    parser = argparse.ArgumentParser()
    parser.add_argument("-dist", type=int, default=0, help = "Min sensor spacing (km)")
    parser.add_argument("-gs", type=int, default=400, help = "Grid size")
    parser.add_argument("-pr", help="Definition of protected region units in meters")
    parser.add_argument("-ap", help = "Definition of antenna patterns unit in Km.")
    parser.add_argument("-anneal", default=False, action="store_true", help="Whether or not to anneal the result")
    parser.add_argument("-of",default="output", help = "Output file name prefix")
    args = parser.parse_args()
    protection_region = args.pr
    do_anneal = args.anneal
    coverage_file = args.ap
    #min_ctr_dist = args.dist
    min_ctr_dist = args.dist
    output_file = args.of
    grid_size = args.gs
    antennacover.NDIVISIONS=grid_size
    
    # Load up the data.
    with open (protection_region, "r") as myfile:
        data=myfile.readlines()

    # Units are in Km.
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
    cover = antennacover.min_antenna_area_cover_greedy(possible_centers, interference_contour, coverage_file, min_center_distance=min_ctr_dist)
    printcover.printAntennaCover(output_file, interference_contour, possible_centers, cover,coverage_file,60,min_ctr_dist)
    if do_anneal:
        annealr = simannealer.SimAnneal(interference_contour, possible_centers, coverage_file,cover)
        annealr.anneal()
        cover = annealr.get_result()
        printcover.printAntennaCover(output_file + "Anneal", interference_contour, possible_centers, cover,coverage_file,60,min_ctr_dist)

    # Return a list of tuples that indicates the coverage.
    sensor_centers = list(set([c[0] for c in cover]))
    sensor_loc_x = [r[0] for r in sensor_centers]
    sensor_loc_y = [r[1] for r in sensor_centers]
    
    f = open(output_file + ".txt","w")
    center_x = [r[0][0] for r in cover]
    center_y = [r[0][1] for r in cover]
    indexes  = [r[1] for r in cover]
    angles  = [r[2] for r in cover]
    f.write(str(center_x) + "\n")
    f.write(str(center_y) + "\n")
    f.write(str(indexes) + "\n")
    f.write(str(angles) + "\n")
    f.write(str(sensor_loc_x) + "\n")
    f.write(str(sensor_loc_y) + "\n")
    f.write(str(len(sensor_centers)) + "\n")
    print "Sensor Count = " , str(len(sensor_centers))
    centers = [c[0] for c in cover]
    sea_excess_area,land_excess_area = excessarea.compute_excess_area_for_antenna_cover(indexes, angles, centers, coverage_file,
            possible_centers, interference_contour)
    f.write(str(sea_excess_area) + "\n")
    f.write(str(land_excess_area) + "\n")
    f.close()

