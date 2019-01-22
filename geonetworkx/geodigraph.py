import networkx as nx
from .geograph import GeoGraph
import geonetworkx as gnx

class GeoDiGraph(GeoGraph, nx.DiGraph):

    def copy(self, as_view=False):
        graph = nx.DiGraph.copy(self, as_view)
        return GeoDiGraph(graph, **self.get_spatial_keys())

    def to_undirected(self, reciprocal=False, as_view=False):
        """Return an undirected copy of the graph."""
        if as_view:
            return nx.DiGraph.to_undirected(self, reciprocal, as_view)
        else:
            graph_class = self.to_undirected_class()
            undirected_graph = nx.DiGraph.to_undirected(self, reciprocal, as_view)
            return graph_class(undirected_graph, **self.get_spatial_keys())

    def to_undirected_class(self):
        """Returns the class to use for empty undirected copies."""
        return gnx.GeoGraph

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph."""
        if as_view:
            return nx.DiGraph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.DiGraph.to_directed(self, as_view)
            return graph_class(directed_graph, **self.get_spatial_keys())

    def to_directed_class(self):
        """Returns the class to use for empty directed copies."""
        return GeoDiGraph


    """ 
    TODO
    def reverse():
        call nx.DiGraph reverse
        Reverse all linestrings
    """