

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

# THIS FILE IS DEPRECATED -- to be removed soon 

import osr,ogr
from shapely.geometry import LineString,Polygon

class Projection :

    def __init__(self,basemap):
        XYSpatialRef = osr.SpatialReference()
        XYSpatialRef.ImportFromEPSG(102008)    # North America Albers Equal Area Conic
        #XYSpatialRef.ImportFromEPSG(3857)    # WGS84 Web Mercator (Auxiliary Sphere)
        # WGS 84 spatial reference
        LatLonSpatialRef = osr.SpatialReference()
        LatLonSpatialRef.ImportFromEPSG(4326)
        # Convert from (longitude, latitude) to projected coordinates
        self.transform = osr.CoordinateTransformation(LatLonSpatialRef, XYSpatialRef)
        # Convert from projected coordinates to (longitude, latitude)
        self.inv_transform = osr.CoordinateTransformation(XYSpatialRef, LatLonSpatialRef)

    def polygon_to_latlon(self,polygon):
        xy_coords = list(polygon.exterior.coords)
        xcoords = [xy_coords[i][0] for i in range(0,len(xy_coords))]
        ycoords = [xy_coords[i][1] for i in range(0,len(xy_coords))]
        zcoords = [0 for i in range(0,len(xy_coords))]
        coast_lons, coast_lats = self.inv_transform.TransformPoints(zip(xcoords,ycoords))
        return LineString(zip(coast_lons,coast_lats,zcoords))

    def lon_lat_to_xy(self,lon,lat):
        xy = self.transform.TransformPoint(lon, lat)
        return xy

    def lons_lats_to_xy(self,lon,lat):
        return self.transform.TransformPoints(zip(lon, lat))
