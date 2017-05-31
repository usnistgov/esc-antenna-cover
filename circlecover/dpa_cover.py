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
	



if __name__=="__main__":
    coast_coords = []
    # East Coast
    #coast_coords.append( [(-80.859375,25.165173),(-66.884766,44.840291)])
    # Gulf Coast
    #coast_coords.append( [(-97.646484,25.958045),(-80.947266,25.005973)])
    # West Coast
    #coast_coords.append( [(-117.070312,32.620870),(-124.892578,48.4000325)])
	
    # Put the lat_0 lon_0 at the geographic center of the USA

    #-125.0011, 24.9493, -66.9326, 49.5904 Centroid:   -95.9669, 37.1669

    # Lambert Conformal Conic Mapping Used for many new USGS maps created after 1957. It replaced the Polyconic projection.
    basemap = Basemap(projection = 'hammer', llcrnrlon = -125.0011, llcrnrlat = 24.9493, urcrnrlon = -66.9326, urcrnrlat = 49.5904, resolution = 'l', lat_0= 37.1669, lon_0=-95.9669)
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="KML file for DPA coordinates")
    parser.add_argument("-f", help="Antenna coverage file")
    args = parser.parse_args()
    dpa_file_name = args.k
    detection_coverage_file = args.f
    if dpa_file_name is None or detection_coverage_file is None:
        parser.print_help()
        sys.exit()
        
    with open(dpa_file_name, 'r') as f:
        datastring = f.read()

    kmlParser = kml.KML()

    kmlDoc = kmlParser.from_string(datastring)
    ns = kmlParser.ns

    coast_coords_xy = pts_coast(basemap)
    counter = 0
    K = 0
    dpa_counter = 1
    lobeId = 0
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
                    dpa_centers.append((x,y))
            elif feature1.name == "East DPAs":
                for feature2 in feature1.features():
                    # east_dpa_10km_1
                    if type(feature2) is fastkml.kml.Document and feature2.name.startswith("east_dpa_10km_19"):
                        testName = "DPA : " + feature2.name
                        print "Processing " + testName
                        for feature3 in feature2.features():
                            latitudes = [p[1] for p in feature3.geometry.coords]
                            longitudes = [p[0] for p in feature3.geometry.coords]
                            x,y = basemap(longitudes, latitudes)
                            dpa_locs = zip(x,y)
                            dpa_polygon = Polygon(dpa_locs)
                            # Determine the DPA center
                            min_dpa_center = min(dpa_centers, key =  lambda t : Point(t).distance(dpa_polygon))
                            for k in range(0,100):
                                # In some cases, 10 KM buffer does not result in any points.
                                # 10 KM buffer do not consider ponts INSIDE the dpa to be candidate locations.
                                extended_dpa_polygon = Polygon(dpa_locs).buffer((1+k*0.5)*10*1000) 
                                candidate_locs = [p for p in coast_coords_xy if extended_dpa_polygon.contains(Point(p)) and not dpa_polygon.contains(Point(p))]
                                if len(candidate_locs) < 20: 
                                    K = K + 1 + k*0.5
                                    counter = counter+1
                                else:
                                    K = K + 1
                                    counter = counter + 1
                                    break

                            print len(candidate_locs)
                            candidate_locs = sorted(candidate_locs,key = lambda t : Point(t).distance(dpa_polygon))

                            if len(candidate_locs) >= 10  :
                                 candidate_locs = candidate_locs[0:10]

                            candidate_locs.append(min_dpa_center)

                        assert len(candidate_locs) > 0
                        #print candidate_locs
                        antennacover.NDIVISIONS = 100
                        cover = antennacover.min_antenna_area_cover_greedy(candidate_locs, dpa_polygon, detection_coverage_file, min_center_distance=0,tol=.001,coverage_units="m")

                        annealer = simannealer.SimAnneal(dpa_polygon,detection_coverage_file,cover,steps = 1000,tol=.001,coverage_units="m")
                        annealer.anneal()
                        newcover = annealer.get_result()
                        #newcover = cover

                        f = kml.Folder(ns, 'Antennas_' + feature2.name, 'Antennas', 'Antenna Placement')
                        feature2.append(f)
                        for c in newcover:
                            center = c[0]
                            index = c[1]
                            angle = c[2]
                            lobe = antennacover.translate_and_rotate(antenna_cover_patterns,center,index,angle)
                            p = kml.Placemark(ns, "antenna"+str(lobeId), 'antenna', 'Index ' + str(index) + " Angle " + str(angle))
                            p.geometry = xy_to_latlon(basemap,lobe)
                            p.name = feature2.name
                            f.append(p)
                            lobeId = lobeId + 1
                        dpa_counter = dpa_counter + 1
                        print "Done " + testName

    print "**** DONE ************"

    with open("output.kml", 'w') as f:
        output = kmlParser.to_string()
        f.write(output)
    
    
    print "average K ", float(K)/float(counter)
                    
                    
            
