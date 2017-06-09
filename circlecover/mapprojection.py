
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

from shapely.geometry import LineString, Polygon

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
# Projected coordinate system

class Projection :

    def __init__(self,basemap):
        self.basemap = basemap

    def lon_lat_to_xy(self,lon,lat):
        x,y = self.basemap(lon,lat)
        return (x,y)
   

    def lons_lats_to_xy(self,lon,lat):
        x,y = self.basemap(lon,lat)
        return zip(x,y)


    def polygon_to_latlon(self,polygon):
        xy_coords = list(polygon.exterior.coords)
        xcoords = [xy_coords[i][0] for i in range(0,len(xy_coords))]
        ycoords = [xy_coords[i][1] for i in range(0,len(xy_coords))]
        zcoords = [0 for i in range(0,len(xy_coords))]
        coast_lons, coast_lats = self.basemap(xcoords,ycoords, inverse = True)
        return LineString(zip(coast_lons,coast_lats,zcoords))

    def xy_to_lonlat(self,x,y):
        lon, lat = self.basemap(x,y, inverse = True)
        return lon,lat
