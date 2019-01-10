import networkx as nx
from .geograph import GeoGraph


class GeoDiGraph(GeoGraph, nx.DiGraph):

    def copy(self, as_view=False):
        graph = nx.DiGraph.copy(self, as_view)
        return GeoDiGraph(graph, **self.get_spatial_keys())

    """
    TODO
    def to_undirected():
        call DiGraph to undirected
        convert to GeoGraph at the end
    """

    """ 
    TODO
    def reverse():
        call nx.DiGraph reverse
        Reverse all linestrings
    """