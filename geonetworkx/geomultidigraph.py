import networkx as nx
from .geomultigraph import GeoMultiGraph
from .geodigraph import GeoDiGraph

class GeoMultiDiGraph(GeoMultiGraph, GeoDiGraph, nx.MultiDiGraph):

    def copy(self, as_view=False):
        graph = nx.MultiDiGraph.copy(self, as_view)
        return GeoMultiDiGraph(graph, **self.get_spatial_keys())

    """
    TODO
    def to_undirected():
        call MultiDiGraph to undirected
        convert to GeoMultiGraph at the end
    """
