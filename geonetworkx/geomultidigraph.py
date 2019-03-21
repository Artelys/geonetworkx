# -*- coding: utf-8 -*-
import networkx as nx
from .geomultigraph import GeoMultiGraph
from .geodigraph import GeoDiGraph
import geonetworkx as gnx


class GeoMultiDiGraph(GeoMultiGraph, GeoDiGraph, nx.MultiDiGraph):
    """
    A directed geographic graph class that can store multiedges.
    """

    def to_nx_class(self):
        """Return the closest networkx class (in the inheritance graph)."""
        return nx.MultiDiGraph

    def to_undirected(self, reciprocal=False, as_view=False):
        """Return an undirected copy of the graph (see ``networkx.MultiDiGraph.to_undirected``)."""
        if as_view:
            return nx.MultiDiGraph.to_undirected(self, reciprocal, as_view)
        else:
            graph_class = self.to_undirected_class()
            undirected_graph = nx.MultiDiGraph.to_undirected(self, reciprocal, as_view)
            return graph_class(undirected_graph)

    def to_undirected_class(self):
        """Returns the class to use for empty undirected copies (see ``networkx.MultiDiGraph.to_undirected_class``)."""
        return gnx.GeoMultiGraph

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph (see ``networkx.MultiDiGraph.to_directed``)."""
        if as_view:
            return nx.MultiDiGraph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.MultiDiGraph.to_directed(self, as_view)
            return graph_class(directed_graph)

    def to_directed_class(self):
        """Returns the class to use for empty directed copies (see ``networkx.MultiDiGraph.to_directed_class``)."""
        return GeoMultiDiGraph
