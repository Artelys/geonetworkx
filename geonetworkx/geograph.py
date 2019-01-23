"""Base class for geographic graphs"""
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point
import geonetworkx as gnx
import geonetworkx.settings as settings


class GeoGraph(nx.Graph):

    def __init__(self, incoming_graph_data=None, **attr):
        self.parse_input_keys(attr)
        super(GeoGraph, self).__init__(incoming_graph_data, **attr)
        self.check_nodes_validity()

    def parse_input_keys(self, attr):
        self.x_key = attr.pop('x_key', settings.X_DEFAULT_KEY)
        self.y_key = attr.pop('y_key', settings.Y_DEFAULT_KEY)
        self.edges_geometry_key = attr.pop('edges_geometry_key', settings.EDGES_GEOMETRY_DEFAULT_KEY)
        self.crs = attr.pop('crs', settings.DEFAULT_CRS)

    def check_nodes_validity(self):
        for n, node_data in self.nodes(data=True):
            if self.x_key not in node_data or self.y_key not in node_data:
                raise ValueError("Unable to find (x, y) coordinates for node: '%s'" % str(n))

    def get_node_coordinates(self, node_name):
        node_data = self.nodes[node_name]
        return [node_data[self.x_key], node_data[self.y_key]]

    def get_nodes_coordinates(self):
        return {n: self.get_node_coordinates(n) for n in self.nodes}

    def get_node_as_point(self, node_name):
        return Point(self.get_node_coordinates(node_name))

    def get_nodes_as_points(self):
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
        return {'x_key': self.x_key,
                'y_key': self.y_key,
                'edges_geometry_key': self.edges_geometry_key,
                'crs': self.crs}

    def set_nodes_coordinates(self, coordinates: dict):
        for n, coords in coordinates.items():
            node_data = self.nodes[n]
            node_data[self.x_key] = coords[0]
            node_data[self.y_key] = coords[1]


    def copy(self, as_view=False):
        graph = nx.Graph.copy(self, as_view)
        return GeoGraph(graph, **self.get_spatial_keys())

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph."""
        if as_view:
            return nx.Graph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.Graph.to_directed(self, as_view)
            return graph_class(directed_graph, **self.get_spatial_keys())

    def to_directed_class(self):
        """Returns the class to use for empty directed copies."""
        return gnx.GeoDiGraph

    def to_undirected(self, as_view=False):
        """Return an undirected copy of the graph."""
        if as_view:
            return nx.Graph.to_undirected(self, as_view)
        else:
            graph_class = self.to_undirected_class()
            undirected_graph = nx.Graph.to_undirected(self, as_view)
            return graph_class(undirected_graph, **self.get_spatial_keys())

    def to_undirected_class(self):
        """Returns the class to use for empty undirected copies."""
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