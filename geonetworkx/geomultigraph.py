# -*- coding: utf-8 -*-
import networkx as nx
import geonetworkx as gnx


class GeoMultiGraph(gnx.GeoGraph, nx.MultiGraph):
    """A undirected geographic graph class that can store multiedges."""

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

    def add_edge(self, u_for_edge, v_for_edge, key=None, **attr):
        """Add a single edge.

        This method exists only for reflecting nx.MultiGraph method so that the multiple inheritance scheme works.

        Examples
        --------
        >>> import geonetworkx as gnx
        >>> g = gnx.GeoMultiGraph()
        >>> g.add_edge(1, 2, 0, geometry=gnx.LineString([(5, 4), (2, 7)]))
        0
        >>> print(g.nodes[1]["geometry"])
        POINT (5 4)
        """
        u_geometry, v_geometry = self._get_nodes_geometries_from_edge_geometry(u_for_edge, v_for_edge,
                                                                               attr.get(self.edges_geometry_key, None))
        result = self.to_nx_class().add_edge(self, u_for_edge, v_for_edge, key, **attr)
        if u_geometry is not None:
            self.nodes[u_for_edge][self.nodes_geometry_key] = u_geometry
        if v_geometry is not None:
            self.nodes[v_for_edge][self.nodes_geometry_key] = v_geometry
        return result
