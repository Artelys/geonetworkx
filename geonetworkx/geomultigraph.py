import networkx as nx
from .geograph import GeoGraph
import geonetworkx as gnx

class GeoMultiGraph(GeoGraph, nx.MultiGraph):

    def copy(self, as_view=False):
        graph = nx.MultiGraph.copy(self, as_view)
        return GeoMultiGraph(graph)

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph."""
        if as_view:
            return nx.MultiGraph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.MultiGraph.to_directed(self, as_view)
            return graph_class(directed_graph, **self.get_spatial_keys())

    def to_directed_class(self):
        """Returns the class to use for empty directed copies."""
        return gnx.GeoMultiDiGraph