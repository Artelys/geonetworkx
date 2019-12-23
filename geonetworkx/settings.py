# -*- coding: utf-8 -*-
# default CRS to set when creating graphs
WGS84_CRS = {'init': 'epsg:4326'}
DEFAULT_CRS = None

# Distance measurement methods. This setting is set in the code.
DISTANCE_MEASUREMENT_METHODS = dict()

# Graph default spatial keys
NODES_GEOMETRY_DEFAULT_KEY = "geometry"
EDGES_GEOMETRY_DEFAULT_KEY = "geometry"


# Intersection prefix for spatial merges
INTERSECTION_PREFIX = "intersect_"
ORIGINAL_EDGE_KEY = "original_edge"
BOUNDARY_NODE_PREFIX = "boundary_"

# read/write convention
EDGE_FIRST_NODE_COLUMN_NAME = "u"
EDGE_SECOND_NODE_COLUMN_NAME = "v"
NODE_ID_COLUMN_NAME = "id"
KNOWN_FILES_EXTENSION = {'GeoJSON': '.geojson', 'GPKG': '.gpkg', 'ESRI Shapefile': '.shp'}
