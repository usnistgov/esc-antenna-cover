import fastkml
import pdb
from fastkml import kml
import sys
import antennacover
import math
import simannealer

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
import argparse



# The desired projection is set with the projection keyword. Default is cyl. Supported values for the projection keyword are:
# Value		Description
# cea		Cylindrical Equal Area
# mbtfpq	McBryde-Thomas Flat-Polar Quartic
# aeqd		Azimuthal Equidistant
# sinu		Sinusoidal
# poly		Polyconic
# omerc		Oblique Mercator
# gnom		Gnomonic
# moll		Mollweide
# lcc		Lambert Conformal
# tmerc		Transverse Mercator
# nplaea	North-Polar Lambert Azimuthal
# gall		Gall Stereographic Cylindrical
# npaeqd	North-Polar Azimuthal Equidistant
# mill		Miller Cylindrical
# merc		Mercator
# stere		Stereographic
# eqdc		Equidistant Conic
# rotpole	Rotated Pole
# cyl		Cylindrical Equidistant
# npstere	North-Polar Stereographic
# spstere	South-Polar Stereographic
# hammer	Hammer
# geos		Geostationary
# nsper		Near-Sided Perspective
# eck4		Eckert IV
# aea		Albers Equal Area
# kav7		Kavrayskiy VII
# spaeqd	South-Polar Azimuthal Equidistant
# ortho		Orthographic
# cass		Cassini-Soldner
# vandg		van der Grinten
# laea		Lambert Azimuthal Equal Area
# splaea	South-Polar Lambert Azimuthal
# robin		Robinson

# For most map projections, the map projection region can either be specified by setting these keywords:

# Keyword	Description
# llcrnrlon	longitude of lower left hand corner of the desired map domain (degrees).
# llcrnrlat	latitude of lower left hand corner of the desired map domain (degrees).
# urcrnrlon	longitude of upper right hand corner of the desired map domain (degrees).
# urcrnrlat	latitude of upper right hand corner of the desired map domain (degrees).
	
# or these

# Keyword	Description
# width		width of desired map domain in projection coordinates (meters).
# height	height of desired map domain in projection coordinates (meters).
# lon_0		center of desired map domain (in degrees).
# lat_0		center of desired map domain (in degrees).

# You can change the resolution of boundary database to use. Can be c (crude), l (low), i (intermediate), h (high), f (full). 
# Coastline data is from the GSHHS (http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html). 
# State, country and river datasets from the Generic Mapping Tools (http://gmt.soest.hawaii.edu).

# Following are the coastlines of the US.

def pts_coast(basemap):
    coast_xy = []
    coastlines = basemap.drawcoastlines()
	# Get the coordinates of coastlines
    coordinates = np.vstack(coastlines.get_segments())
    coast_lons, coast_lats = basemap(coordinates[:,0],coordinates[:,1], inverse = True)
    x, y = basemap(coast_lons, coast_lats)
    coast_xy = coast_xy + zip(x,y)
    return coast_xy


def xy_to_latlon(basemap,polygon):
    xy_coords = list(polygon.exterior.coords)
    xcoords = [xy_coords[i][0] for i in range(0,len(xy_coords))]
    ycoords = [xy_coords[i][1] for i in range(0,len(xy_coords))]
    zcoords = [0 for i in range(0,len(xy_coords))]
    coast_lons, coast_lats = basemap(xcoords,ycoords, inverse = True)
    return LineString(zip(coast_lons,coast_lats,zcoords))



def distance(p1,p2):
    return math.sqrt(math.pow((float(p1[0]) - float(p2[0])),2) + math.pow((float(p1[1]) - float(p2[1])),2))
	


def parse_forbidden_region(basemap,forbidden_file_name):
    with open(forbidden_file_name, 'r') as f:
        datastring = f.read()

    poly = None
    kmlParser = kml.KML()
    kmlDoc = kmlParser.from_string(datastring)
    for feature in kmlParser.features():
        for feature1 in feature.features():
            poly = feature1.geometry
            coords = poly.exterior.coords
            longs = [coords[i][0] for i in range(0,len(coords))]
            lats =  [coords[i][1] for i in range(0,len(coords))] 
            x,y = basemap(longs,lats, inverse = False)
            return Polygon(zip(x,y))
            
            

def is_forbidden_loc(loc,forbidden_polygons):
    for polygon in forbidden_polygons:
        if polygon.contains(Point(loc)):
            print "istVerboten "
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
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="KML file for DPA coordinates")
    parser.add_argument("-f", help="Antenna coverage file")
    parser.add_argument("-e", default=None, help="Sensor placement exclusion region")
    args = parser.parse_args()
    dpa_file_name = args.k
    detection_coverage_file = args.f
    forbidden_region_files = args.e
    # Regions where we MAY NOT place sensors:
    forbidden_polygons = []
    if forbidden_region_files is not None:
        forbidden_region_file_names = forbidden_region_files.split(",")
        for i in range(0,len(forbidden_region_file_names)):
            poly = parse_forbidden_region(basemap,forbidden_region_file_names[i])
            forbidden_polygons.append(poly)
        
        
    if dpa_file_name is None or detection_coverage_file is None:
        parser.print_help()
        sys.exit()
        
    with open(dpa_file_name, 'r') as f:
        datastring = f.read()

    kmlParser = kml.KML()

    kmlDoc = kmlParser.from_string(datastring)

    coast_coords_xy = pts_coast(basemap)
    counter = 0
    K = 0
    dpa_counter = 1
    lobeId = 0
    sensorId = 0
    antenna_cover_patterns = antennacover.read_detection_coverage(detection_coverage_file,coverage_units="m")
    folders = []
    dpa_centers = []
    for feature in kmlParser.features():
        for feature1 in feature.features():
            print "feature1.name = " + feature1.name
            if feature1.name == "East DPA Centers" :
               for feature2 in feature1.features():
                    print feature2.name
                    lon = feature2.geometry.x
                    lat = feature2.geometry.y
                    x,y = basemap(lon,lat)
                    dpa_centers.append((feature2.name,(x,y)))
            elif feature1.name == "East DPAs":
                for feature2 in feature1.features():
                    # east_dpa_10km_1
                    if type(feature2) is fastkml.kml.Document and feature2.name.startswith("east_dpa_10km"):
                        testName = "DPA : " + feature2.name
                        print "Processing " + testName
                        for feature3 in feature2.features():
                            latitudes = [p[1] for p in feature3.geometry.coords]
                            longitudes = [p[0] for p in feature3.geometry.coords]
                            x,y = basemap(longitudes, latitudes)
                            dpa_locs = zip(x,y)
                            dpa_polygon = Polygon(dpa_locs)
                            # Determine the DPA center
                            min_dpa_name,min_dpa_center = min(dpa_centers, key =  lambda t : Point(t[1]).distance(dpa_polygon))
                            for k in range(0,100):
                                # In some cases, 10 KM buffer does not result in any points or in very few points.
                                # we keep extending the boundary till we get 20 points to choose from where we can place sensors.
                                extended_dpa_polygon = Polygon(dpa_locs).buffer((1+k*0.5)*10*1000) 
                                candidate_locs = [p for p in coast_coords_xy if extended_dpa_polygon.contains(Point(p)) and not dpa_polygon.contains(Point(p))]
                                if len(candidate_locs) < 5: 
                                    K = K + 1 + k*0.5
                                    counter = counter+1
                                else:
                                    K = K + 1
                                    counter = counter + 1
                                    break

                            candidate_locs = [ cand for cand in candidate_locs if not is_forbidden_loc(cand,forbidden_polygons) ]
                            

                            # Pick the top 10 candidate locations. 
                            candidate_locs = sorted(candidate_locs,key = lambda t : Point(t).distance(dpa_polygon))

                            if len(candidate_locs) >= 5  :
                                 candidate_locs = candidate_locs[0:4]

                            # Add our dpa center to the candidate locations (not sure why).
                            candidate_locs.append(min_dpa_center)

                        assert len(candidate_locs) > 0
                        #print candidate_locs
                        antennacover.NDIVISIONS = 100
                        cover = antennacover.min_antenna_area_cover_greedy(candidate_locs, dpa_polygon, detection_coverage_file, min_center_distance=0,tol=.001,coverage_units="m")

                        annealer = simannealer.SimAnneal(dpa_polygon,detection_coverage_file,cover,steps = 1000,tol=.001,coverage_units="m")
                        annealer.anneal()
                        newcover = annealer.get_result()
                        #newcover = cover

                        f = kml.Folder(kmlParser.ns, 'Antennas_' + min_dpa_name, 'Antennas', 'Antenna Angles')
                        feature2.append(f)
                        for c in newcover:
                            center = c[0]
                            index = c[1]
                            angle = c[2]
                            lobe = antennacover.translate_and_rotate(antenna_cover_patterns,center,index,angle)
                            p = kml.Placemark(kmlParser.ns, "antenna"+str(lobeId), 'antenna', 'Index ' + str(index) + " Angle " + str(angle))
                            p.geometry = xy_to_latlon(basemap,lobe)
                            p.name = feature2.name
                            f.append(p)
                            lobeId = lobeId + 1
                        
                        f = kml.Folder(kmlParser.ns, 'Sensors_' + min_dpa_name, 'Sensors', 'Antenna Placement')
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
                            p = kml.Placemark(kmlParser.ns,"sensor_" + str(sensorId), "sensor:" + min_dpa_name, "Index " + str(indexes[i]))
                            p.geometry = Point(lon,lat)
                            p.styleUrl = "#msn_shaded_dot1"
                            i = i + 1
                            sensorId = sensorId + 1
                            f.append(p)

                        dpa_counter = dpa_counter + 1
                        print "Done " + testName

    print "**** DONE ************"

    with open("output.kml", 'w') as f:
        output = kmlParser.to_string(prettyprint=True)
        f.write(output)
    
    
    print "average K ", float(K)/float(counter)
                    
                    
            
