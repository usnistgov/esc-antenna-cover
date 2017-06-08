import fastkml
import pdb
from fastkml import kml
import sys
import antennacover
import math
import simannealer
import excessarea
import os
import re

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import MultiPoint
import argparse
from mapprojection import Projection
#from epsgprojection import Projection



def pts_coast(projection,basemap):
    coast_xy = []
    coastlines = basemap.drawcoastlines()
	# Get the coordinates of coastlines
    coordinates = np.vstack(coastlines.get_segments())
    coast_lons, coast_lats = basemap(coordinates[:,0],coordinates[:,1], inverse = True)
    vals = projection.lons_lats_to_xy(coast_lons, coast_lats)
    coast_xy = coast_xy + vals
    return coast_xy


def distance(p1,p2):
    return math.sqrt(math.pow((float(p1[0]) - float(p2[0])),2) + math.pow((float(p1[1]) - float(p2[1])),2))
	


def parse_forbidden_region(porjection,forbidden_file_name):
    with open(forbidden_file_name, 'r') as f:
        datastring = f.read()

    poly = None
    kmlDoc = kml.KML()
    kmlDoc.from_string(datastring)
    for feature in kmlDoc.features():
        for feature1 in feature.features():
            poly = feature1.geometry
            coords = poly.exterior.coords
            lons = [coords[i][0] for i in range(0,len(coords))]
            lats =  [coords[i][1] for i in range(0,len(coords))] 
            return Polygon(projection.lons_lats_to_xy(lons,lats))
            
            

def is_forbidden_loc(loc,forbidden_polygons):
    for polygon in forbidden_polygons:
        if polygon.contains(Point(loc)):
            print "placement ist verboten!"
            return True
    return False


if __name__=="__main__":
    
    # The following is not used explicitly in the code but is here for
    # informational purposes.
    # boundaries of the US East Coast
    # [(-80.859375,25.165173),(-66.884766,44.840291)]
    # boundaries of the US Gulf Coast
    # [(-97.646484,25.958045),(-80.947266,25.005973)]
    # boundaries of the US West Coast
    # [(-117.070312,32.620870),(-124.892578,48.4000325)]

    # Put the lat_0 lon_0 at the geographic center of the USA. Here are the extents of USA
    #-125.0011, 24.9493, -66.9326, 49.5904 Centroid:   -95.9669, 37.1669

    # Lambert Conformal Conic Mapping Used for many new USGS maps created after 1957. 
    # It replaced the Polyconic projection. However, we want to use
    # an area preserving projection "hammer" and kav7 work well to preserve areas.

    basemap = Basemap(projection = 'hammer', llcrnrlon = -125.0011, llcrnrlat = 24.9493, urcrnrlon = -66.9326, urcrnrlat = 49.5904, resolution = 'l', lat_0= 37.1669, lon_0=-95.9669)

    projection = Projection(basemap)

    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="KML file for DPA coordinates")
    parser.add_argument("-f", help="Antenna coverage files")
    parser.add_argument("-e", default=None, help="Sensor placement exclusion region")
    parser.add_argument("-d", default = None, help = "DPA id regular expression for which to compute cover (eg. east_dpa_10km_* for all east_dpa)")
    args = parser.parse_args()
    dpa_file_name = args.k
    dpa_file_path = os.path.dirname(os.path.abspath(dpa_file_name))
    detection_coverage_dir = args.f
    forbidden_region_files = args.e
    if args.d is not None:
        dpa_name = re.compile(args.d)
    else:
        dpa_name = None
    # Regions where we MAY NOT place sensors:
    forbidden_polygons = []
    if forbidden_region_files is not None:
        forbidden_region_file_names = forbidden_region_files.split(",")
        for i in range(0,len(forbidden_region_file_names)):
            poly = parse_forbidden_region(projection,forbidden_region_file_names[i])
            forbidden_polygons.append(poly)
        
        
    if dpa_file_name is None or detection_coverage_dir is None:
        parser.print_help()
        sys.exit()
        
    with open(dpa_file_name, 'r') as f:
        datastring = f.read()

    kmlDoc = kml.KML()

    kmlDoc.from_string(datastring)

    coast_coords_xy = pts_coast(projection,basemap)
    counter = 0
    K = 0
    sensor_counter = 0
    dpa_counter = 1
    lobeId = 0
    sensorId = 0
    folders = []
    dpa_centers = []
    dpa_covers = []
    dpa_polygons = []
    total_area = 0
    for feature in kmlDoc.features():
        for feature1 in feature.features():
            print "feature1.name = " + feature1.name
            if feature1.name == "East DPA Centers" or feature1.name == "West DPA Centers":
               for feature2 in feature1.features():
                    print feature2.name
                    lon = feature2.geometry.x
                    lat = feature2.geometry.y
                    xy = projection.lon_lat_to_xy(lon,lat)
                    dpa_centers.append((feature2.name,xy))
            elif feature1.name == "East DPAs" or feature1.name == "West DPAs":
                for feature2 in feature1.features():
                    # east_dpa_10km_1
                    if type(feature2) is fastkml.kml.Document and dpa_name is None or dpa_name.match(feature2.name):
                        testName = "DPA : " + feature2.name

                        sensorPlacement = kml.KML()
                        ns = '{http://www.opengis.net/kml/2.2}'

                        # Create a KML Document and add it to the KML root object
                        placementDoc = kml.Document(ns, 'sensor_placement', 'Sensor Placement', 'Sensor Placement for DPA : ' + feature2.name)
                        sensorPlacement.append(placementDoc)

                        print "Processing " + testName
                        for feature3 in feature2.features():
                            latitudes = [p[1] for p in feature3.geometry.coords]
                            longitudes = [p[0] for p in feature3.geometry.coords]
                            x,y = basemap(longitudes, latitudes)
                            dpa_locs = zip(x,y)
                            # BUGBUG -- some polygons in the DPA are self intersecting. (e.g. dpa 29)
                            if not Polygon(dpa_locs).is_simple:
                                dpa_polygon = Polygon(MultiPoint(dpa_locs).convex_hull)
                            else:
                                dpa_polygon = Polygon(dpa_locs)
                            dpa_locs = dpa_polygon.exterior.coords
                            for k in range(0,500):
                                # In some cases, 10 KM buffer does not result in any points or in very few points.
                                # we keep extending the boundary till we get 20 points to choose from where we can place sensors.
                                extended_dpa_polygon = Polygon(dpa_locs).buffer((1+k*0.5)*10*1000) 
                                candidate_locs = [p for p in coast_coords_xy if extended_dpa_polygon.contains(Point(p)) and not dpa_polygon.contains(Point(p))]
                                if len(candidate_locs) < 100: 
                                    K = K + 1 + k*0.5
                                    counter = counter+1
                                else:
                                    K = K + 1
                                    counter = counter + 1
                                    break

                            # Exclude candidate locations that are in the forbidden regions
                            candidate_locs = [ cand for cand in candidate_locs if not is_forbidden_loc(cand,forbidden_polygons) ]
                            
                            # Sort by distance from the closest DPA.
                            candidate_locs = sorted(candidate_locs,key = lambda t : Point(t).distance(dpa_polygon))

                            # Pick the top 20 candidate locations that are closest to the DPA 
                            if len(candidate_locs) >= 30  :
                                 candidate_locs = candidate_locs[0:29]



                        assert len(candidate_locs) > 0
                        antennacover.NDIVISIONS = 100

                        if dpa_polygon.area == 0:
                            print "Empty dpa"
                            continue

                        if os.path.isdir(detection_coverage_dir) :
                            # Form a tuple list of aperture angle and detection coverage file name.
                            # We want to apply increasing angles if there is more than one sensor in a DPA.
                            detection_coverage_file_list = [(antennacover.read_aperture_angle(detection_coverage_dir+ "/" + dcf),detection_coverage_dir + "/" + dcf) 
                                                for dcf in os.listdir(detection_coverage_dir)]

                            # Sort these in order of increasing antenna angle.
                            sorted_detection_coverage_files = sorted(detection_coverage_file_list, key = lambda(t) : t[0])
                            # Put the result in a list
                            detection_coverage_files = [s[1]  for s in sorted_detection_coverage_files]
                        else:
                            # We are given a single antenna cover file.
                            detection_coverage_files = [detection_coverage_dir]
                
                            
                        print "Antenna Detection Coverage files ", detection_coverage_files
                        for detection_coverage_file in detection_coverage_files:
                            aperture_angle = antennacover.read_aperture_angle(detection_coverage_file)
                            antenna_cover_patterns = antennacover.read_detection_coverage(detection_coverage_file,coverage_units="m")
                            cover = antennacover.min_antenna_area_cover_greedy(candidate_locs, dpa_polygon, detection_coverage_file, min_center_distance=0,tol=.005,coverage_units="m")

                            if len(cover) > 1:
                                # There isn't any point in annealing if you have only one lobe.
                                annealer = simannealer.SimAnneal(dpa_polygon,detection_coverage_file,cover,steps = 1000,tol=.005,coverage_units="m")
                                annealer.anneal()
                                newcover = annealer.get_result()
                                if len(newcover) == 1:
                                    break
                            else:
                                newcover = cover
                                break

                        # Keep the cover around for later analysis.
                        dpa_covers.append(newcover)
                        dpa_polygons.append(dpa_polygon)
                        f = kml.Folder(kmlDoc.ns, 'Antennas_' + feature2.name, 'Antennas', 'Antenna Angles')
                        placementDoc.append(f)
                        feature2.append(f)
                        for c in newcover:
                            center = c[0]
                            index = c[1]
                            angle = c[2]
                            lobe = antennacover.translate_and_rotate(antenna_cover_patterns,center,index,angle)
                            angle_degrees = (angle/math.pi*180)%360
                            p = kml.Placemark(kmlDoc.ns, "antenna"+str(lobeId), 'antenna', 'Sensitivity (dBm): ' + str(antenna_cover_patterns[index].sensitivity_dbm) 
                                            + "\naperture angle " + str(aperture_angle)
                                            + "\nAngle: " + str(angle)  + " rad " + str(angle_degrees) + " degrees. " )
                            p.geometry = projection.xy_to_latlon(lobe)
                            p.name = feature2.name
                            f.append(p)
                            lobeId = lobeId + 1

                        interference_contour = list(dpa_polygon.exterior.coords)
                        indexes = [newcover[i][1] for i in range(0,len(newcover))]
                        angles  = [newcover[i][2] for i in range(0,len(newcover))]
                        cover_centers = [newcover[i][0] for i in range(0,len(newcover))]
                        
                        sea_excess_area,land_excess_area,outage_area,coverage_area = excessarea.compute_excess_area_for_antenna_cover(indexes, 
                                        angles, cover_centers, detection_coverage_file, candidate_locs, interference_contour,units="m")

                        
                        f = kml.Folder(kmlDoc.ns, 'Sensors_' + feature2.name, 'Sensors', 'Antenna Placement')
                        placementDoc.append(f)
                        feature2.append(f)
                        sensor_locs = list(set([c[0] for c in newcover]))
                        indexes = []
                        for sensor_loc in sensor_locs:
                            for c in newcover:
                                if c[0] == sensor_loc:
                                    indexes.append(c[1])
                                    break

                        for i in range(0,len(sensor_locs)) :
                            sensor_loc = sensor_locs[i]
                            lon,lat = basemap(sensor_loc[0],sensor_loc[1],inverse=True)
                            p = kml.Placemark(kmlDoc.ns,"sensor_" + str(sensorId), "sensor: " + feature2.name, 
                                            "Sensitivity (dBm): " + str(antenna_cover_patterns[indexes[i]].sensitivity_dbm) + "\nlon : " + str(lon) + " lat : " +  str(lat))
                            p.geometry = Point(lon,lat)
                            p.styleUrl = "#msn_shaded_dot1"
                            i = i + 1
                            sensorId = sensorId + 1
                            f.append(p)

                        placementDoc.description = placementDoc.description + \
                                                    "\nDPA Name " + feature2.name +\
                                                    "\nsensor_count " + str(len(sensor_locs)) + \
                                                    "\nsea_excess_area " + str(sea_excess_area) + \
                                                    "\nland_excess_area " + str(land_excess_area) + \
                                                    "\noutage_area " + str(outage_area) + \
                                                    "\ntotal_coverage_area " + str(coverage_area)

                        dpa_counter = dpa_counter + 1
                        total_area = total_area + dpa_polygon.area
                        sensor_counter = sensor_counter + len(sensor_locs)
                        
                        with open(dpa_file_path + "/" + feature2.name + ".kml", 'w') as f:
                            output = sensorPlacement.to_string(prettyprint=True)
                            f.write(output)
                        print "Done " + testName + " Sensor_count " + str(len(sensor_locs))


    with open(dpa_file_path + "/antenna-cover.kml","w") as f:
        output = kmlDoc.to_string(prettyprint=True)
        f.write(output)
        
    
    print "Computing intersection area "
    intersection_area = 0
    for i in range(0,len(dpa_covers)):
        for j in range(0,len(dpa_covers)):
            if j != i :
                # The composite polygon is the union of all the lobes.
                centers = [c[0] for c in dpa_covers[i]   ]
                indexes = [c[1] for c in dpa_covers[i]   ]
                angles =  [c[2] for c in dpa_covers[i]  ]
                cover_polygons = excessarea.generate_antenna_cover_polygons(indexes,angles,centers,antenna_cover_patterns)
                union = cover_polygons[0]
                for k in range(1,len(cover_polygons)):
                    union = union.union(cover_polygons[k])
                diff = dpa_polygons[j].difference(dpa_polygons[j].difference(union)).area
                intersection_area = intersection_area + diff

    ratio = intersection_area / total_area
    print "intersection_area " , intersection_area, " total_area ", total_area
    print "Ratio of intersection to total area " , str(ratio)

    with open("analysis.txt",'w') as f:
        f.write("intersection_area " + str(intersection_area) +  " total_area " + str(total_area) + "\n")
        f.write("Ratio of intersection to total area "  +  str(ratio) + " (probability of 2 DPA's being activated with a single sensor)\n")
        f.write("Sensor count " + str(sensor_counter) + "\n")
            
                
    print "**** DONE ************"

    
    
                    
                    
            
