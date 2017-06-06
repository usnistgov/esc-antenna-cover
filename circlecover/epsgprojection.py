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

    def xy_to_latlon(self,polygon):
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
