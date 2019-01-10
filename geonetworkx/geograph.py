"""Base class for geographic graphs"""
import networkx as nx

class GeoGraph(nx.Graph):
    EDGES_GEOMETRY_DEFAULT_KEY = "geometry"
    X_DEFAULT_KEY = 'x'
    Y_DEFAULT_KEY = 'y'

    def __init__(self, incoming_graph_data=None, **attr):
        self.parse_input_keys(**attr)
        super(GeoGraph, self).__init__(incoming_graph_data, **attr)
        self.check_nodes_validity()

    def parse_input_keys(self, **attr):
        if 'x_key' in attr:
            self.x_key = attr['x_key']
            del attr['x_key']
        else:
            self.x_key = GeoGraph.X_DEFAULT_KEY
        if 'y_key' in attr:
            self.y_key = attr['y_key']
            del attr['y_key']
        else:
            self.y_key = GeoGraph.Y_DEFAULT_KEY
        if 'edges_geometry_key' in attr:
            self.edges_geometry_key = attr['edges_geometry_key']
            del attr['edges_geometry_key']
        else:
            self.edges_geometry_key = GeoGraph.EDGES_GEOMETRY_DEFAULT_KEY

    def check_nodes_validity(self):
        for n, node_data in self.nodes(data=True):
            if self.x_key not in node_data or self.y_key not in node_data:
                raise ValueError("Unable to find (x, y) coordinates for node: '%s'" % str(n))

    def get_node_coordinates(self, node_name):
        node_data = self.nodes[node_name]
        return [node_data[self.x_key], node_data[self.y_key]]

    def get_spatial_keys(self):
        return {'x_key': self.x_key,
                'y_key': self.y_key,
                'edges_geometry_key': self.edges_geometry_key}

    def copy(self, as_view=False):
        graph = nx.Graph.copy(self, as_view)
        return GeoGraph(graph, **self.get_spatial_keys())
    """
    TODO
    def to_directed():
        call nx.Graph to_directed
        convert to GeoDiGraph at the end
    """

