import networkx as nx
from .geograph import GeoGraph


class GeoMultiGraph(GeoGraph, nx.MultiGraph):

    def copy(self, as_view=False):
        graph = nx.MultiGraph.copy(self, as_view)
        return GeoMultiGraph(graph)
    """
    TODO
    def to_directed():
        call nx.MultiGraph to_directed
        convert to GeoMultiDiGraph at the end
    """