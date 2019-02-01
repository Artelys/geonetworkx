"""Base class for geographic graphs"""
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point, LineString
import geonetworkx as gnx
import geonetworkx.settings as settings


class GeoGraph(nx.Graph):
    """
    This class extends the `networkx.Graph` to represent a graph that have a geographical meaning. Nodes are located
    with their coordinates (x, y) and edges can be represented with a given broken line (using
    `shapely.geometry.LineString` objects). Each graph has its own keys for naming nodes x and y coordinates and edges
    geometry (`x_key`, `y_key`, `edges_geometry_key`). A coordinate reference system (CRS) can be defined for a graph
    and will be used for some methods managing earth coordinates (especially for distances). For now, the only supported
    CRS is the WGS84 standard (EPSG:4326). All nodes must have defined coordinates, otherwise an error will be raised.
    """
    def __init__(self, incoming_graph_data=None, **attr):
        self.x_key = attr.pop('x_key', settings.X_DEFAULT_KEY)
        self.y_key = attr.pop('y_key', settings.Y_DEFAULT_KEY)
        self.edges_geometry_key = attr.pop('edges_geometry_key', settings.EDGES_GEOMETRY_DEFAULT_KEY)
        self.crs = attr.pop('crs', settings.DEFAULT_CRS)
        super(GeoGraph, self).__init__(incoming_graph_data, **attr)
        self.check_nodes_validity()

    def check_nodes_validity(self):
        """Check that all nodes have x and y coordinates."""
        for n, node_data in self.nodes(data=True):
            if self.x_key not in node_data or self.y_key not in node_data:
                raise ValueError("Unable to find (x, y) coordinates for node: '%s'" % str(n))

    def get_node_coordinates(self, node_name):
        """Return the coordinates of a given node."""
        node_data = self.nodes[node_name]
        return [node_data[self.x_key], node_data[self.y_key]]

    def get_nodes_coordinates(self):
        """Return all nodes coordinates within a dictionary."""
        return {n: self.get_node_coordinates(n) for n in self.nodes}

    def get_node_as_point(self, node_name):
        """Return a node as a `shapely.geometry.Point` object."""
        return Point(self.get_node_coordinates(node_name))

    def get_nodes_as_points(self):
        """Return all nodes as `shapely.geometry.Point` objects within a dictionary."""
        return {n: self.get_node_as_point(n) for n in self.nodes}

    def get_nodes_as_point_series(self):
        """Return the nodes as a `geopandas.GeoSeries` of `shapely.geometry.Point`."""
        nodes_as_points = self.get_nodes_as_points()
        point_series = gpd.GeoSeries(nodes_as_points)
        point_series.crs = self.crs
        return point_series

    def get_edges_as_line_series(self):
        """Return the edges as a `geopandas.GeoSeries` of `shapely.geometry.LineString`."""
        lines = nx.get_edge_attributes(self, self.edges_geometry_key)
        line_series = gpd.GeoSeries(lines)
        line_series.crs = self.crs
        return line_series

    def get_spatial_keys(self):
        """Return the current graph spatial keys."""
        return {'x_key': self.x_key,
                'y_key': self.y_key,
                'edges_geometry_key': self.edges_geometry_key,
                'crs': self.crs}

    def set_nodes_coordinates(self, coordinates: dict):
        """Set nodes coordinates with a given dictionary of coordinates (can be used for a subset of all nodes)."""
        for n, coords in coordinates.items():
            node_data = self.nodes[n]
            node_data[self.x_key] = coords[0]
            node_data[self.y_key] = coords[1]

    def to_nx_class(self):
        """Return the closest networkx class (in the inheritance graph)."""
        return nx.Graph

    def copy(self, as_view=False):
        """Return a copy of the graph (see `networkx.Graph.copy`)."""
        nx_graph_class = self.to_nx_class()
        graph = nx_graph_class.copy(self, as_view)
        return self.__class__(graph, **self.get_spatial_keys())

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph (see `networkx.Graph.to_directed`)."""
        if as_view:
            return nx.Graph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.Graph.to_directed(self, as_view)
            return graph_class(directed_graph, **self.get_spatial_keys())

    def to_directed_class(self):
        """Returns the class to use for empty directed copies (see `networkx.Graph.to_directed_class`)."""
        return gnx.GeoDiGraph

    def to_undirected(self, as_view=False):
        """Return an undirected copy of the graph (see `networkx.Graph.to_undirected`)."""
        if as_view:
            return nx.Graph.to_undirected(self, as_view)
        else:
            graph_class = self.to_undirected_class()
            undirected_graph = nx.Graph.to_undirected(self, as_view)
            return graph_class(undirected_graph, **self.get_spatial_keys())

    def to_undirected_class(self):
        """Returns the class to use for empty undirected copies (see `networkx.Graph.to_undirected_class`)."""
        return gnx.GeoGraph

    def to_crs(self, crs=None, epsg=None, inplace=False):
        """Transform edge geometries and nodes coordinates to a new coordinate reference system."""
        if self.crs is None:
            raise ValueError('Cannot transform naive geometries. Please set a crs on the graph first.')
        if inplace:
            graph = self
        else:
            graph = self.copy()
        # Get new nodes coordinates
        nodes_as_points = graph.get_nodes_as_point_series()
        transformed_nodes = nodes_as_points.to_crs(crs, epsg)
        # Get new edges coordinates
        edges_as_lines = graph.get_edges_as_line_series()
        transformed_edges = edges_as_lines.to_crs(crs, epsg)
        # Operate the transformation
        for n, point in transformed_nodes.iteritems():
            node_data = graph.nodes[n]
            node_data[graph.x_key] = point.x
            node_data[graph.y_key] = point.y
        for e, line in transformed_edges.iteritems():
            edge_data = graph.edges[e]
            edge_data[graph.edges_geometry_key] = line
        graph.crs = transformed_nodes.crs
        if not inplace:
            return graph


    def nodes_to_gdf(self) -> gpd.GeoDataFrame:
        """
        Create a `geopandas.GeoDataFrame` from nodes of the current graph. The 'geometry' attribute is used for shapes
        writing (from `geopandas.GeoDataFrame.DEFAULT_GEO_COLUMN_NAME`).

        :return: The resulting GeoDataFrame : one row is a node
        """
        nodes = {node: data for node, data in self.nodes(data=True)}
        gdf_nodes = gpd.GeoDataFrame(nodes).T
        gdf_nodes[settings.NODE_ID_COLUMN_NAME] = gdf_nodes.index
        gdf_nodes[settings.GPD_GEOMETRY_KEY] = gdf_nodes.apply(lambda row: Point(row[self.x_key], row[self.y_key]),
                                                               axis=1)
        gdf_nodes.crs = self.crs
        return gdf_nodes

    def edges_to_gdf(self) -> gpd.GeoDataFrame:
        """
        Create a `geopandas.GeoDataFrame` from edges of the current graph. The 'geometry' attribute is used for shapes
        writing (from `geopandas.GeoDataFrame.DEFAULT_GEO_COLUMN_NAME`).

        :return: The resulting GeoDataFrame : one row is an edge
        """
        # create a list to hold our edges, then loop through each edge in the graph
        edges = []
        for u, v, data in self.edges(data=True):
            # for each edge, add key and all attributes in data dict to the edge_details
            edge_details = {settings.EDGE_FIRST_NODE_COLUMN_NAME: u, settings.EDGE_SECOND_NODE_COLUMN_NAME: v}
            for attr_key in data:
                edge_details[attr_key] = data[attr_key]
            # if edge doesn't already have a geometry attribute, create one now
            if self.edges_geometry_key not in data:
                point_u = Point((self.nodes[u][self.x_key], self.nodes[u][self.y_key]))
                point_v = Point((self.nodes[v][self.x_key], self.nodes[v][self.y_key]))
                edge_details[settings.GPD_GEOMETRY_KEY] = LineString([point_u, point_v])
            else:
                line = edge_details[self.edges_geometry_key]
                del edge_details[self.edges_geometry_key]
                edge_details[settings.GPD_GEOMETRY_KEY] = line
            edges.append(edge_details)
        # create a GeoDataFrame from the list of edges and set the CRS
        gdf_edges = gpd.GeoDataFrame(edges)
        gdf_edges.crs = self.crs
        return gdf_edges
