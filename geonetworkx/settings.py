# default CRS to set when creating graphs
WGS84_CRS = {'init': 'epsg:4326'}
USED_CRS = WGS84_CRS
DEFAULT_CRS = None

# default discretization parameter for extended projection methods
DISCRETIZATION_TOLERANCE = 1e-4

# Graph default spatial keys
NODES_GEOMETRY_DEFAULT_KEY = "geometry"
EDGES_GEOMETRY_DEFAULT_KEY = "geometry"


# Intersection prefix for spatial merges
INTERSECTION_PREFIX = "intersect_"
ORIGINAL_EDGE_KEY = "original_edge"

# read/write convention
EDGE_FIRST_NODE_COLUMN_NAME = "u"
EDGE_SECOND_NODE_COLUMN_NAME = "v"
NODE_ID_COLUMN_NAME = "id"
KNOWN_FILES_EXTENSION = {'GeoJSON': '.geojson', 'GPKG': '.gpkg', 'ESRI Shapefile': '.shp'}