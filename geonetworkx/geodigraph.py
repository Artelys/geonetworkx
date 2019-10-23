# -*- coding: utf-8 -*-
import networkx as nx
from .geograph import GeoGraph
import geonetworkx as gnx


class GeoDiGraph(GeoGraph, nx.DiGraph):
    """Base class for directed geographic graphs.
    
    Because edges are directed, it supposes that the edges lines are well-ordered. Namely, that the first point of the
    line matches with the coordinates of the first vertex of the edge (or is at least close) and vice versa with the
    last point of the line and the second. If this is not the case, the method ``order_well_lines`` can be useful to
    make sure of that.
    """

    def to_nx_class(self):
        """ """
        return nx.DiGraph

    def to_undirected(self, reciprocal=False, as_view=False):
        """Return an undirected copy of the graph (see ``networkx.DiGraph.to_undirected``)."""
        if as_view:
            return nx.DiGraph.to_undirected(self, reciprocal, as_view)
        else:
            graph_class = self.to_undirected_class()
            undirected_graph = nx.DiGraph.to_undirected(self, reciprocal, as_view)
            return graph_class(undirected_graph)

    def to_undirected_class(self):
        """Returns the class to use for empty undirected copies (see ``networkx.DiGraph.to_undirected_class``)."""
        return gnx.GeoGraph

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph (see ``networkx.DiGraph.to_directed``)."""
        if as_view:
            return nx.DiGraph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.DiGraph.to_directed(self, as_view)
            return graph_class(directed_graph)

    def to_directed_class(self):
        """Returns the class to use for empty directed copies (see ``networkx.DiGraph.to_directed_class``)."""
        return GeoDiGraph

    """
    TODO
    def reverse():
        call nx.DiGraph reverse
        Reverse all linestrings
    """
