# default CRS to set when creating graphs
WGS84_CRS = {'init': 'epsg:4326'}
USED_CRS = WGS84_CRS
DEFAULT_CRS = None

# default discretization parameter for extended projection methods
DISCRETIZATION_TOLERANCE = 1e-4

# Geopandas "geometry" attribute
GPD_GEOMETRY_KEY = "geometry"

# Graph default spatial keys
EDGES_GEOMETRY_DEFAULT_KEY = "geometry"
X_DEFAULT_KEY = 'x'
Y_DEFAULT_KEY = 'y'


# Intersection prefix for spatial merges
INTERSECTION_PREFIX = "intersect_"
