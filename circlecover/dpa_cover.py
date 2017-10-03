
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
import string

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import MultiPoint
from shapely.ops import polygonize, unary_union
import argparse
from mapprojection import Projection
from pykml import parser as kml_parser
#from epsgprojection import Projection

MAX_CANDIDATE_LOCS=30
DISTANCE_FROM_COAST = 50000 # distance from coast to consider a point inside sea in meters

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


    # Lambert Conformal Conic Mapping Used for many new USGS maps created after 1957. 
    # It replaced the Polyconic projection. However, we want to use
    # Unfortunately Conic Mappings map to a cone and therefore do not preserve areas.
    # an area preserving projection "hammer" and kav7 work well to preserve areas.
    # Put the lat_0 lon_0 at the geographic center of the USA. Here are the extents of USA
    #-125.0011, 24.9493, -66.9326, 49.5904 Centroid:   -95.9669, 37.1669

    basemap = Basemap(projection = 'hammer', llcrnrlon = -125.0011, llcrnrlat = 24.9493, urcrnrlon = -66.9326, urcrnrlat = 49.5904, resolution = 'l', lat_0= 37.1669, lon_0=-95.9669)

    projection = Projection(basemap)

    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="KML file for DPA coordinates")
    parser.add_argument("-f", help="Antenna coverage file or directory containing multiple cover files.\n" + 
                                    "If a directory is specified then antenna covers are tried in  \n " + 
                                    "increasing order of  Aperture until a single sensor can cover the \n" + 
                                    "entire DPA OR the largest aperture has been reached.")
    parser.add_argument("-e", default=None, help="Optionally specified sensor placement forbidden regions (comma separated).")
    parser.add_argument("-d", default = "east_dpa*", help = "DPA id regular expression for which to compute cover " +
                                                    "(eg. east_dpa* for all east_dpa).")
    parser.add_argument("-o",default="./", help="Output file path to directory where you want the resultant kml " + 
                                                "files to be written (default is current directory)")
    args = parser.parse_args()
    dpa_file_name = args.k
    dpa_file_path = os.path.dirname(os.path.abspath(dpa_file_name))
    detection_coverage_dir = args.f
    forbidden_region_files = args.e
    output_directory = args.o
    
    # East DPA and West DPA are best computed independently.
    if not args.d.startswith("east_dpa") and not args.d.startswith("west_dpa") :
        print "DPA identifier should start with east or west"
        parser.print_help()
        sys.exit()

    print "---- args.d = ", args.d
    dpa_regexp = re.compile(args.d)

    # Make a directory to keep the results.
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Regions where we MAY NOT place sensors:
    forbidden_polygons = []
    if forbidden_region_files is not None:
        forbidden_region_file_names = forbidden_region_files.split(",")
        for i in range(0,len(forbidden_region_file_names)):
            poly = parse_forbidden_region(projection,forbidden_region_file_names[i])
            forbidden_polygons.append(poly)
        
        
    if dpa_file_name is None or detection_coverage_dir is None:
        print "Detection coverage file name or directory must be specified "
        parser.print_help()
        sys.exit()
        
    # Read the DPA definition file.
    with open(dpa_file_name, 'r') as f:
        datastring = f.read()

    kmlDoc = kml.KML()

    # Parse the DPA definition file.
    kmlDoc.from_string(datastring)
    feature = list(kmlDoc.features())
    print "---- features ", len(feature), feature[0].name
    feature1 = list(feature[0].features())
    print "---- features1.name ", len(feature1), feature1[0].name

    # Get the coastal coordinates from the basemap
    coast_coords_xy = pts_coast(projection,basemap)
    coast_coords_linestring = LineString(coast_coords_xy)
    # Keeps track of the number of sensors.
    sensor_counter = 0
    # Keep track of the number of antenna lobes.
    lobe_counter = 0
    folders = []
    # The dinstinguished markers for each DPA ( these are not used)
    dpa_centers = []
    # The covers representing the DPAs.
    dpa_covers = []
    # the polygons representing the DPAs.
    dpa_polygons = []
    # Each element of this array is a union of polygons that covers the corresponding dpa_polygon
    antenna_lobes = []
    total_area = 0
    
    root = kml_parser.fromstring(datastring)
    print "---- all", root.Folder.name.text
    #dpa_id = root.Folder.Document.name.text
    #print "dpa_id ", root.Folder.Document.name.text
    feature_index = 0
    dpa_metrics =[]
    for Doc in root.Folder.Document:
        #print "---- coord ", Doc.Placemark.Polygon.outerBoundaryIs.LinearRing.coordinates
        dpa_id = Doc.name.text
        print "---- dpa_id ", dpa_id
        if dpa_regexp.match(dpa_id):
            #print "------- matched dpa ------ ", dpa_id
            testName = "DPA : " + dpa_id
            print "Processing " + dpa_id
            sensorPlacement = kml.KML()
            ns = '{http://www.opengis.net/kml/2.2}'

            # Create a KML Document and add it to the KML root object
            placementDoc = kml.Document(ns, 'sensor_placement', 'Sensor Placement', 'Sensor Placement for DPA : ' + dpa_id)
            sensorPlacement.append(placementDoc)
            p = Doc.Placemark.Polygon.outerBoundaryIs.LinearRing.coordinates.text
            coord = string.split(p)
            #print "------ coord = ", coord
            
            latitudes =[]
            longitudes =[]
            for p in string.split(Doc.Placemark.Polygon.outerBoundaryIs.LinearRing.coordinates.text):
                pts = string.split(p, ',')
                longitudes.append(float(pts[0]))
                latitudes.append(float(pts[1]))
            #print "---- lats  ", latitudes
            #latitudes = [string.split(p, ',') for p in string.split(Doc.Placemark.Polygon.outerBoundaryIs.LinearRing.coordinates.text)]
            #longitudes = [p[0] for p in Doc.Placemark.Polygon.outerBoundaryIs.LinearRing.coordinates]
            x,y = basemap(longitudes, latitudes)
            dpa_locs = zip(x,y)
            #some polygons in the DPA are self intersecting. (e.g. dpa 29) 
            #so we split these into polygons and take the largest polygon
            #pdb.set_trace()
            ls = LineString(dpa_locs)     
            mls = unary_union(ls)
            max_area = -1
            dpa_polygon = None
            # Use the biggest polygon 
            # Polygonize returns an iterator of polygons formed from the lines.
            # Pick the biggest polygon and use it to represent the DPA
            for p in polygonize(mls):
                if p.area > max_area:
                    dpa_polygon = p.buffer(0)
                    max_area = p.area
            
            dist = []
            print "len(coast_coord) = ", len(coast_coords_xy)
            for p in coast_coords_xy:
                # Exclude candidate locations that are in the forbidden regions
                if is_forbidden_loc(p,forbidden_polygons):
                    dist.append(4294967295) # make the distance large
                    continue
                point = Point(p)
                #print "---- point -------" , point
                dist.append(dpa_polygon.distance(point))
            
            dist = sorted(range(len(dist)), key=lambda k: dist[k]) # dist has the indices of coast_coord_xy sorted in the order of distance from the DPA
            candidate_locs = []
            for i in range(MAX_CANDIDATE_LOCS):
                candidate_locs.append(coast_coords_xy[dist[i]])

            ## ----- start of code to print candidate location in long lat (this code is not needed)
            lon = []
            lat = []
            for i in range(len(candidate_locs)):
               p = Point(candidate_locs[i])
               cand_lon, cand_lat = basemap(p.x, p.y, inverse=True)
               lon.append(cand_lon)
               lat.append(cand_lat)

            lonlat = zip(lon, lat)
            #print "lonlat of candidate loc ", lonlat
            ## ----- end of code to print candidate location in long lat

            print "------ polygon area ------ ", dpa_polygon.area
            """
            dpa_locs = dpa_polygon.exterior.coords
            for k in range(0,500):
                # In some cases, 10 KM buffer does not result in any points or in very few points.
                # we keep extending the boundary till we get 20 points to choose from where we can place sensors.
                extended_dpa_polygon = dpa_polygon.buffer((1+k*0.5)*1*1000) 
                candidate_locs = [p for p in coast_coords_xy if extended_dpa_polygon.contains(Point(p)) and not dpa_polygon.contains(Point(p))]
                if len(candidate_locs) >= 100: 
                    break

            # Exclude candidate locations that are in the forbidden regions
            candidate_locs = [ cand for cand in candidate_locs if not is_forbidden_loc(cand,forbidden_polygons) ]
            
            # Sort by distance from the closest DPA.
            candidate_locs = sorted(candidate_locs,key = lambda t : Point(t).distance(dpa_polygon))

            # HACK ALERT Pick the top 30 candidate locations that are closest to the DPA in terms of distance
            if len(candidate_locs) >= 30  :
                candidate_locs = candidate_locs[0:29]
            """


            assert len(candidate_locs) > 0
            #print "---- candidate locs ----", candidate_locs

            if dpa_polygon.area == 0:
                print "Empty dpa"
                continue

            # Set global parameter (this is ugly -- need to change this)
            antennacover.NDIVISIONS = 100
        
            # Read the detection coverage files. First check if we are given a directory with multiple
            # detection coverage files - we'll try these one at a time.
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
            min_cover_lobes = None
            for detection_coverage_file in detection_coverage_files:
                # Read the aperture angle from the file.
                aperture_angle = antennacover.read_aperture_angle(detection_coverage_file)
                antenna_cover_patterns = antennacover.read_detection_coverage(detection_coverage_file,coverage_units="m")
                cover = antennacover.min_antenna_area_cover_greedy(candidate_locs, dpa_polygon, detection_coverage_file, min_center_distance=0,tol=.005,coverage_units="m")
                cover_lobes = None
                for c in cover:
                    center = c[0]
                    index = c[1]
                    angle = c[2]
                    lobe = antennacover.translate_and_rotate(antenna_cover_patterns,center,index,angle)
                    if cover_lobes is None:
                        cover_lobes = lobe
                    else:
                        cover_lobes = cover_lobes.union(lobe)
                #print "--- AREA ---", cover_lobes.area
                if  (min_cover_lobes is None):
                    min_cover_lobes = cover_lobes
                    min_cover = cover
                elif (cover_lobes.area < min_cover_lobes.area):
                    min_cover_lobes = cover_lobes
                    min_cover = cover

            #print "--- COVER AREA ---", min_cover_lobes.area
            cover = min_cover
            cover_lobes = min_cover_lobes
            if len(cover) > 1:
               # There isn't any point in annealing if you have only one lobe.
               annealer = simannealer.SimAnneal(dpa_polygon,detection_coverage_file,cover,steps = 1000,tol=.005,coverage_units="m")
               annealer.anneal()
               newcover = annealer.get_result()
               if len(newcover) == 1:
                   print "newcover has one element"
            else:
                print "----- newcover has one element"
                newcover = cover


            # ****** TAKE THIS OUT WHEN YOU PUT SIM ANNEALING ****
            #newcover = cover
            # Keep the cover around for later analysis. dpa_covers is the cover for each DPA
            dpa_covers.append(newcover)
            # dpa_polygons is an array containing DPA polygons
            dpa_polygons.append(dpa_polygon)
            # Update the KML document.
            f = kml.Folder(ns, 'Antennas_' + dpa_id, 'Antennas', 'Antenna Angles')
            placementDoc.append(f)
            feature1[feature_index].append(f)
            for c in newcover:
                center = c[0]
                index = c[1]
                angle = c[2]
                lobe = antennacover.translate_and_rotate(antenna_cover_patterns,center,index,angle)
                # Note that our frame of reference is Due EAST but the conventional 
                # way of specifying angles is due north.
                angle_degrees = (angle/math.pi*180 -90)%360 
                p = kml.Placemark(ns, "antenna"+str(lobe_counter), 'antenna', 'Sensitivity (dBm): ' + str(antenna_cover_patterns[index].sensitivity_dbm) 
                                  + "\nAperture angle: " + str(aperture_angle) + " degrees" 
                                  + "\nAzimuth angle: " +  str(float(np.round(angle_degrees,2))) + " degrees from due North. " )
                p.geometry = projection.polygon_to_latlon(lobe)
                p.name = Doc.name.text
                f.append(p)
                lobe_counter = lobe_counter + 1

            print "---- covered lobes area, lobe_counter ----", cover_lobes.area, lobe_counter
            # antenna_lobes is an array consisting of composite antenna cover polygons.
            antenna_lobes.append(cover_lobes)
            interference_contour = list(dpa_polygon.exterior.coords)
            indexes = [newcover[i][1] for i in range(0,len(newcover))]
            angles  = [newcover[i][2] for i in range(0,len(newcover))]
            cover_centers = [newcover[i][0] for i in range(0,len(newcover))]

            # Compute the excess area. We cant use shapely for this directly because the
            # union could be non-simple. 
            # Get the bounds of the area of interest.
            minx,miny,maxx,maxy = cover_lobes.union(dpa_polygon).bounds
            outage_counter = 0
            excess_area_counter = 0
            coverage_counter = 0
            sea_excess_counter = 0
            # Try to keep the integration regions square in shape.
            aspect_ratio = (maxy-miny)/(maxx - minx)
            if aspect_ratio < 1:
                YDIVS = 100
                XDIVS = int(100/aspect_ratio)
            else:
                XDIVS = 100
                YDIVS = int(100/aspect_ratio)
            
            integration_element_area = float((maxx-minx)*(maxy-miny) / ( XDIVS * YDIVS ))
            for i in range(0,XDIVS):
                for j in range(0,YDIVS):
                    x = minx + i*(maxx-minx)/XDIVS
                    y = miny + j*(maxy-miny)/YDIVS
                    p = Point(x,y)
                    if dpa_polygon.contains(p) and not cover_lobes.contains(p):
                        outage_counter = outage_counter + 1
                    if cover_lobes.contains(p) and not dpa_polygon.contains(p):
                        excess_area_counter = excess_area_counter + 1 
                        # excess area outside the DPA; includes area on the land as well
                        if (not basemap.is_land(x,y)):
                            sea_excess_counter += 1 # this includes area from nbring DPAs
                         # p.distance(coast_coords_linestring) does not give correct results; so avoid using it.
#                        if ((not basemap.is_land(x,y)) and (p.distance(coast_coords_linestring) > DISTANCE_FROM_COAST)) :
                    if dpa_polygon.contains(p):
                        coverage_counter = coverage_counter + 1


            # Compute the areas for THIS DPA. 
            outage_area = outage_counter*integration_element_area
            # Note that the Excess area here does not distinguish between sea and land.
            # Also does not distinguish between this DPA and its neighbor. 
            # We evaluate that later across all DPAs.
            excess_area = excess_area_counter*integration_element_area
            coverage_area = coverage_counter*integration_element_area
            sea_excess_area = sea_excess_counter*integration_element_area
            # includes area covered in the nbring DPAs
            print "----- ", outage_counter, excess_area_counter, coverage_counter, sea_excess_counter, integration_element_area

            f = kml.Folder(ns, 'Sensors_' + dpa_id, 'Sensors', 'Antenna Placement')
            placementDoc.append(f)
            feature1[feature_index].append(f)
            sensor_locs = list(set([c[0] for c in newcover]))
            indexes = []
            for sensor_loc in sensor_locs:
                for c in newcover:
                    if c[0] == sensor_loc:
                        indexes.append(c[1])
                        break

            # Make up a placemark for each sensor.
            for i in range(0,len(sensor_locs)) :
                sensor_loc = sensor_locs[i]
                lon,lat = basemap(sensor_loc[0],sensor_loc[1],inverse=True)
                p = kml.Placemark(ns,"sensor_" + str(sensor_counter), "sensor: " + dpa_id, 
                                  "Sensitivity (dBm): " + str(antenna_cover_patterns[indexes[i]].sensitivity_dbm) + "\nlon : " + str(lon) + " lat : " +  str(lat))
                p.geometry = Point(lon,lat)
                p.styleUrl = "#msn_shaded_dot1"
                i = i + 1
                # Update the sensor counter.
                sensor_counter = sensor_counter + 1
                f.append(p)

            print "--------- num sensors ------ ", len(sensor_locs)
            placementDoc.description = placementDoc.description + \
                                       "\nDPA Name " + dpa_id +\
                                       "\nsensor_count " + str(len(sensor_locs)) + \
                                       "\nexcess_area (sq. m): " + str(float(np.round(excess_area,2))) + \
                                       "\noutage_area (sq. m): " + str(float(np.round(outage_area,2))) + \
                                       "\ndpa_area (sq. m): " + str(float(np.round(coverage_area,2)))

            total_area = total_area + dpa_polygon.area
            individual_dpa_metrics = {
#                'grid_area' : integration_element_area,
                'dpa_id' : dpa_id,
                'excess_area' : excess_area,
                'outage_area' : outage_area,
                'coverage_area' : coverage_area,
                'dpa_area' : dpa_polygon.area,
                'sea_excess_area' : sea_excess_area, # we need to take out nbring DPA area (later)
                'nbr_excess_area' : 0.0,  # excess area over ndring DPAs; will be calculated later
                'p_outage' : 0.0,
                'p_fa_sea' : 0.0,
                'p_fa_nbrDPA' : 0.0
                }
            dpa_metrics.append(individual_dpa_metrics)
        
            with open(output_directory + "/" + dpa_id + ".kml", 'w') as f:
                output = sensorPlacement.to_string(prettyprint=True)
                f.write(output)
            print "Done " + testName + " Sensor_count " + str(len(sensor_locs))
        feature_index = feature_index + 1    


    with open(output_directory + "/antenna-cover.kml","w") as f:
        output = kmlDoc.to_string(prettyprint=True)
        f.write(output)
        
    
    print "Computing intersection area "
    # Intersecting area is the area outside the DPA that is covered by the sensors of this DPA.
    # It represents the area which is incorrectly detected by a DPA sensors meant for 
    # a given DPA.
    intersecting_area = 0
    for i in range(0,len(dpa_covers)):
        # antenna_lobes[i] is an array consisting of the union of lobes covering DPA i.
        antenna_cover = antenna_lobes[i]
        # note that dpa_coveres, antenna_cover have the same indices 'i'. This is the antenna cover
        # for dpa_polygons[i].
        nbr_excess_area = 0
        for j in range(0,len(dpa_polygons)):
            if j != i :
                # intersecting area between the cover and DPA is the area within some nbring DPA that is covered
                # by a sensor not belonging to the DPA.
                # We compute the intersection of the cover with each dpa 
                # polygon except the one it is supposed to cover.
                common_poly = antenna_cover.intersection(dpa_polygons[j])
                area = common_poly.area
                intersecting_area += area
                nbr_excess_area += area

        dpa_metrics[i]['nbr_excess_area'] = nbr_excess_area
        dpa_metrics[i]['sea_excess_area'] -= nbr_excess_area
        dpa_metrics[i]['p_outage'] = dpa_metrics[i]['outage_area']/dpa_metrics[j]['dpa_area']
        dpa_metrics[i]['p_fa_sea'] = dpa_metrics[i]['sea_excess_area']/dpa_metrics[j]['coverage_area']
        dpa_metrics[i]['p_fa_nbrDPA'] = dpa_metrics[i]['nbr_excess_area']/dpa_metrics[j]['coverage_area']
        print dpa_metrics[i]
        print "----------------------------------------------------------------"

    # Probability of false DPA activation is the probability of activating a DPA not belonging to a given DPA cover.
    """prob_false_dpa_activation = intersecting_area / total_area

    print "intersection_area " , intersecting_area, " total_area ", total_area
    print "Ratio of intersection to total area " , str(prob_false_dpa_activation)

    # Now compute the excess coverage that is not in a DPA that is in the sea. This is 
    # The false alarm ratio. 

    cover_union = antenna_lobes[0]

    for i in range(1, len(antenna_lobes)):
        cover_union = cover_union.union(antenna_lobes[i])

    dpa_polygons_union = dpa_polygons[0]
    
    for i in range(1,len(dpa_polygons)):
        dpa_polygons_union = dpa_polygons_union.union(dpa_polygons[i])

    minx,miny,maxx,maxy = cover_union.union(dpa_polygons_union).bounds
    
    aspect_ratio = (maxy - miny) / (maxx - minx)
    
    if aspect_ratio < 1 :
        NDIVISIONS_Y = 200
        NDIVISIONS_X = int(200/aspect_ratio)
    else:
        NDIVISIONS_Y = int(200/aspect_ratio)
        NDIVISIONS_X =  200

    deltax = (maxx - minx)/NDIVISIONS_X
    deltay =  (maxy - miny)/NDIVISIONS_Y
    sea_excess_area_counter = 0
    total_area_counter = 0
    outage_area_counter = 0

    grid_point_area = abs((maxx-minx)*(maxy-miny)) / (NDIVISIONS_X*NDIVISIONS_Y)
    
    # Get a line string of coastal coordinates.
    for i in range(0,NDIVISIONS_X):
        for j in range(0,NDIVISIONS_Y):
            x,y = minx + i*deltax, miny + j*deltay
            p = Point(x,y)
            # If neither the cover nor the DPA polygon contains the point P then ignore it.
            if cover_union.contains(p) or dpa_polygons_union.contains(p) :
                total_area_counter = total_area_counter + 1
                if not dpa_polygons_union.contains(p):
                    # p is in the cover but NOT in a DPA polygon.
                    # Check if the point is on the sea and beyond 10 km from shore
                    if not basemap.is_land(x,y) and p.distance(coast_coords_linestring) > 10*1000 :
                        sea_excess_area_counter = sea_excess_area_counter + 1
                elif not cover_union.contains(p):
                    # p is in DPA but not in the cover union. This is an outage.
                    outage_area_counter = outage_area_counter + 1


    sea_excess_area = sea_excess_area_counter*grid_point_area
    total_area =  total_area_counter*grid_point_area
    outage_area = outage_area_counter*grid_point_area
    false_detection_probability = np.round(sea_excess_area/total_area,2)
    outage_probability = np.round(outage_area/total_area,2)"""

    

    with open( output_directory + "/analysis.txt",'w') as f:
        f.write("t : total_dpa_area (sq. m) = " + str(float(np.round(total_area,2))) + "\n")
        #f.write("s : Sensor area outside DPA of sensor covered by DPA (total) (sq. m) =  " + str(float(np.round(intersecting_area,2))) + "\n") 
        #f.write("k : Sea Excess Area = " + str(float(np.round(sea_excess_area,2))) + "\n")
        #f.write("h : outage area  = " + str(float(np.round(outage_area,2))) + "\n\n")
        #f.write("Probability redundant DPA activation (s/t)  =  "  +  str(float(np.round(prob_false_dpa_activation,2))) + "\n\n")
        #f.write("Probability of false detection ( sea excess Area to total DPA area (k/t) ) =  " + str(false_detection_probability) + "\n\n")
        #f.write("Probability of missed detection ( outage area to total DPA area (h/t) ) = " + str(outage_probability) + "\n\n")
        f.write("Sensor count =  " + str(sensor_counter) + "\n\n")
        f.write("Antenna count = " + str(lobe_counter) + "\n\n")            

    print "**** DONE ************"

    
    
                    
                    
            
