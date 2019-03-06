"""Base class for geographic graphs"""
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point, LineString
import geonetworkx as gnx
import geonetworkx.settings as settings


class GeoGraph(nx.Graph):
    """
    This class extends the ``networkx.Graph`` to represent a graph that have a geographical meaning. Nodes are located
    with their coordinates (x, y) (using ``shapely.geometry.Point`` objects) and edges can be represented with a given
    broken line (using ``shapely.geometry.LineString`` objects). Each graph has its own keys for naming nodes and edges
    geometry (``nodes_geometry_key``, ``edges_geometry_key``). A coordinate reference system (CRS) can be defined for a
    graph and will be used for some methods managing earth coordinates (especially for distances). For now, the only
    supported CRS is the WGS84 standard (EPSG:4326). All nodes must have defined coordinates, otherwise an error will be
    raised.
    """
    def __init__(self, incoming_graph_data=None, **attr):
        super(GeoGraph, self).__init__(incoming_graph_data, **attr)
        self.check_nodes_validity()

    def check_nodes_validity(self):
        """Check that all nodes have x and y coordinates."""
        for n, node_data in self.nodes(data=True):
            if self.nodes_geometry_key not in node_data:
                raise ValueError("Unable to find geometry for node: '%s'" % str(n))

    @property
    def nodes_geometry_key(self):
        """Attribute name for the edges geometry attributes. This graph attribute appears in the attribute dict G.graph
         keyed by the string `"edges_geometry_key"` as well as an attribute `G.edges_geometry_key`"""
        return self.graph.get('nodes_geometry_key', settings.NODES_GEOMETRY_DEFAULT_KEY)

    @nodes_geometry_key.setter
    def nodes_geometry_key(self, s):
        self.graph['nodes_geometry_key'] = s

    @property
    def edges_geometry_key(self):
        """Attribute name for the edges geometry attributes. This graph attribute appears in the attribute dict G.graph
         keyed by the string `"edges_geometry_key"` as well as an attribute `G.edges_geometry_key`"""
        return self.graph.get('edges_geometry_key', settings.EDGES_GEOMETRY_DEFAULT_KEY)

    @edges_geometry_key.setter
    def edges_geometry_key(self, s):
        self.graph['edges_geometry_key'] = s

    @property
    def crs(self):
        """Coordinate Reference System of the graph. This graph attribute appears in the attribute dict G.graph keyed
        by the string `"crs"` as well as an attribute `G.crs`"""
        return self.graph.get('crs', gnx.DEFAULT_CRS)

    @crs.setter
    def crs(self, c):
        self.graph['crs'] = c

    def get_node_coordinates(self, node_name):
        """Return the coordinates of a given node."""
        point = self.get_node_as_point(node_name)
        return [point.x, point.y]

    def get_nodes_coordinates(self):
        """Return all nodes coordinates within a dictionary."""
        return {n: self.get_node_coordinates(n) for n in self.nodes}

    def get_node_as_point(self, node_name):
        """Return a node as a `shapely.geometry.Point` object."""
        node_data = self.nodes[node_name]
        return node_data[self.nodes_geometry_key]

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
        return {'nodes_geometry_key': self.nodes_geometry_key,
                'edges_geometry_key': self.edges_geometry_key,
                'crs': self.crs}

    def set_nodes_coordinates(self, coordinates: dict):
        """Set nodes coordinates with a given dictionary of coordinates (can be used for a subset of all nodes)."""
        for n, coords in coordinates.items():
            node_data = self.nodes[n]
            node_data[self.nodes_geometry_key] = Point(coords)

    def to_nx_class(self):
        """Return the closest networkx class (in the inheritance graph)."""
        return nx.Graph

    def copy(self, as_view=False):
        """Return a copy of the graph (see `networkx.Graph.copy`)."""
        nx_graph_class = self.to_nx_class()
        graph = nx_graph_class.copy(self, as_view)
        return self.__class__(graph)

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph (see `networkx.Graph.to_directed`)."""
        if as_view:
            return nx.Graph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.Graph.to_directed(self, as_view)
            return graph_class(directed_graph)

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
            return graph_class(undirected_graph)

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
            node_data[graph.nodes_geometry_key] = point
        for e, line in transformed_edges.iteritems():
            edge_data = graph.edges[e]
            edge_data[graph.edges_geometry_key] = line
        graph.crs = transformed_nodes.crs
        if not inplace:
            return graph

    def nodes_to_gdf(self) -> gpd.GeoDataFrame:
        """
        Create a ``geopandas.GeoDataFrame`` from nodes of the current graph. The column representing the geometry is
        named after the current ``nodes_geometry_key`` attribute.

        :return: The resulting GeoDataFrame : one row is a node
        """
        nodes = {node: data for node, data in self.nodes(data=True)}
        gdf_nodes = gpd.GeoDataFrame(nodes).T
        gdf_nodes[settings.NODE_ID_COLUMN_NAME] = gdf_nodes.index
        gdf_nodes.set_geometry(self.nodes_geometry_key, inplace=True)
        gdf_nodes.crs = self.crs
        return gdf_nodes

    def edges_to_gdf(self) -> gpd.GeoDataFrame:
        """
        Create a ``geopandas.GeoDataFrame`` from edges of the current graph. The column representing the geometry is
        named after the current ``edges_geometry_key`` attribute.

        :return: The resulting GeoDataFrame : one row is an edge
        """
        # create a list to hold our edges, then loop through each edge in the graph
        edges = []
        for u, v, data in self.edges(data=True):
            # for each edge, add key and all attributes in data dict to the edge_details
            edge_details = {settings.EDGE_FIRST_NODE_COLUMN_NAME: u, settings.EDGE_SECOND_NODE_COLUMN_NAME: v}
            edge_details.update(data)
            # if edge doesn't already have a geometry attribute, create one now
            if self.edges_geometry_key not in data:
                point_u = self.nodes[u][self.nodes_geometry_key]
                point_v = self.nodes[v][self.nodes_geometry_key]
                edge_details[self.edges_geometry_key] = LineString([point_u, point_v])
            edges.append(edge_details)
        # create a GeoDataFrame from the list of edges and set the CRS
        gdf_edges = gpd.GeoDataFrame(edges)
        gdf_edges.set_geometry(self.edges_geometry_key, inplace=True)
        gdf_edges.crs = self.crs
        return gdf_edges
