# -*- coding: utf-8 -*-
import networkx as nx
from .geograph import GeoGraph
import geonetworkx as gnx


class GeoMultiGraph(GeoGraph, nx.MultiGraph):
    """
    A undirected geographic graph class that can store multiedges.
    """

    def to_nx_class(self):
        """Return the closest networkx class (in the inheritance graph)."""
        return nx.MultiGraph

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph (see ``networkx.MultiGraph.to_directed``)."""
        if as_view:
            return nx.MultiGraph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.MultiGraph.to_directed(self, as_view)
            return graph_class(directed_graph)

    def to_directed_class(self):
        """Returns the class to use for empty directed copies (see ``networkx.MultiGraph.to_directed_class``)."""
        return gnx.GeoMultiDiGraph

    def to_undirected(self, as_view=False):
        """Return an undirected copy of the graph (see ``networkx.MultiGraph.to_undirected``)."""
        if as_view:
            return nx.MultiGraph.to_undirected(self, as_view)
        else:
            graph_class = self.to_undirected_class()
            undirected_graph = nx.MultiGraph.to_undirected(self, as_view)
            return graph_class(undirected_graph)

    def to_undirected_class(self):
        """Returns the class to use for empty undirected copies (see ``networkx.MultiGraph.to_undirected_class``).."""
        return GeoMultiGraph
