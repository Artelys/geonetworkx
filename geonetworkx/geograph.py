# -*- coding: utf-8 -*-
"""Base class for geographic graphs"""
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiPoint
import geonetworkx as gnx
import geonetworkx.settings as settings
import itertools


class GeoGraph(nx.Graph):
    """This class extends the ``networkx.Graph`` to represent a graph that have a geographical meaning. Nodes are
    located with their coordinates (x, y) (using ``shapely.geometry.Point`` objects) and edges can be represented with a
    given broken line (using ``shapely.geometry.LineString`` objects). Each graph has its own keys for naming nodes and
    edges geometry (``nodes_geometry_key``, ``edges_geometry_key``). A coordinate reference system (CRS) can be defined
    for a graph and will be used for some methods managing earth coordinates (especially for distances).
    All nodes must have defined coordinates, otherwise a default coordinates are used.

    Raises
    ------
    ValueError
        If the all nodes don't have valid coordinates.

    See Also
    --------
    networkx.Graph
    GeoDiGraph
    GeoMultiGraph
    GeoMultiDiGraph

    """
    DEFAULT_NODE_GEOMETRY = Point(0, 0)

    def get_default_node_dict(self):
        """Return the default node attribute dictionary."""
        return {self.nodes_geometry_key: self.default_node_geometry}

    def __init__(self, incoming_graph_data=None, **attr):
        self.default_node_geometry = attr.get("default_node_geometry",
                                              self.DEFAULT_NODE_GEOMETRY)
        self.node_attr_dict_factory = self.get_default_node_dict
        super(GeoGraph, self).__init__(incoming_graph_data, **attr)
        self.check_nodes_validity()

    def check_nodes_validity(self):
        """Check that all nodes have geometries."""
        for n, node_data in self.nodes(data=True):
            self.node_attr_dict_check(node_data)

    def node_attr_dict_check(self, attr):
        """Check that the given attribute dictionary contains mandatory fields for a node."""
        if self.nodes_geometry_key not in attr:
            raise ValueError("Node geometry must be in node attributes.")

    @property
    def nodes_geometry_key(self):
        """Attribute name for the edges geometry attributes. This graph attribute appears in the attribute dict
        `G.graph` keyed by the string ``"edges_geometry_key"`` as well as an attribute ``G.nodes_geometry_key``"""
        return self.graph.get('nodes_geometry_key', settings.NODES_GEOMETRY_DEFAULT_KEY)

    @nodes_geometry_key.setter
    def nodes_geometry_key(self, s: str):
        """Sets the node geometry key with the given value"""
        self.graph['nodes_geometry_key'] = s

    @property
    def edges_geometry_key(self):
        """Attribute name for the edges geometry attributes. This graph attribute appears in the attribute dict
        `G.graph` keyed by the string ``"edges_geometry_key"`` as well as an attribute ``G.edges_geometry_key``"""
        return self.graph.get('edges_geometry_key', settings.EDGES_GEOMETRY_DEFAULT_KEY)

    @edges_geometry_key.setter
    def edges_geometry_key(self, s: str):
        """Sets the edges geometry key with the given value"""
        self.graph['edges_geometry_key'] = s

    @property
    def crs(self):
        """Coordinate Reference System of the graph. This graph attribute appears in the attribute dict `G.graph` keyed
        by the string ``"crs"`` as well as an attribute ``G.crs``"""
        return self.graph.get('crs', gnx.DEFAULT_CRS)

    @crs.setter
    def crs(self, c):
        """Sets the current crs"""
        self.graph['crs'] = c

    def get_node_coordinates(self, node_name) -> list:
        """Return the coordinates of the given node.

        Parameters
        ----------
        node_name
            Name of the node on which the coordinates are browsed.

        Returns
        -------
        list
            A two-element list containing (x,y) coordinates of the given node.

        See Also
        --------
        get_nodes_coordinates, get_node_as_point, get_nodes_as_points

        """
        point = self.get_node_as_point(node_name)
        return [point.x, point.y]

    def get_nodes_coordinates(self) -> dict:
        """Return all nodes coordinates within a dictionary.
        
        Returns
        -------
        dict
            Dictionary containing the coordinates of the each node of the graph.

        See Also
        --------
        get_node_coordinates, get_node_as_point, get_nodes_as_points

        """
        return {n: self.get_node_coordinates(n) for n in self.nodes}

    def get_node_as_point(self, node_name):
        """Return a node as a ``shapely.geometry.Point`` object.

        Parameters
        ----------
        node_name :
            Name of the node on which the geometry is browsed.

        Returns
        -------
        shapely.geometry.Point
            The point representing the located node.

        See Also
        --------
        get_node_coordinates, get_node_coordinates, get_nodes_as_points

        """
        node_data = self.nodes[node_name]
        return node_data[self.nodes_geometry_key]

    def get_nodes_as_points(self) -> dict:
        """Return all nodes as ``shapely.geometry.Point`` objects within a dictionary.

        Returns
        -------
        dict
            Dictionary containing the geometry of each node of the graph.

        See Also
        --------
        get_node_coordinates, get_node_coordinates, get_node_as_point

        """
        return {n: self.get_node_as_point(n) for n in self.nodes}

    def get_nodes_as_point_series(self) -> gpd.GeoSeries:
        """Return the nodes as a ``geopandas.GeoSeries`` of ``shapely.geometry.Point``.

        Returns
        -------
        gpd.GeoSeries
            Series containing all nodes geometries. Its CRS is the graph CRS.

        See Also
        --------
        nodes_to_gdf, get_edges_as_line_series

        """
        nodes_as_points = self.get_nodes_as_points()
        point_series = gpd.GeoSeries(nodes_as_points)
        point_series.crs = self.crs
        return point_series

    def get_nodes_as_multipoint(self) -> MultiPoint:
        """Return nodes geometries as a ``shapely.geometry.MultiPoint``.

        Returns
        -------
        MultiPoint
            MutltiPoint containing all nodes geometries.
        """
        return MultiPoint([self.get_node_as_point(n) for n in self.nodes])

    def get_edges_as_line_series(self) -> gpd.GeoSeries:
        """Return the edges as a ``geopandas.GeoSeries`` of ``shapely.geometry.LineString``.

        Returns
        -------
        gpd.GeoSeries
            Series containing all edges geometries. Its CRS is the graph CRS.

        See Also
        --------
        edges_to_gdf, get_nodes_as_point_series

        """
        lines = nx.get_edge_attributes(self, self.edges_geometry_key)
        line_series = gpd.GeoSeries(lines)
        line_series.crs = self.crs
        return line_series

    def get_spatial_keys(self) -> dict:
        """Return the current graph spatial keys.

        Returns
        -------
        dict
            Dictionary containing spatial keys (nodes and edges geometry keys and crs).
        """
        return {'nodes_geometry_key': self.nodes_geometry_key,
                'edges_geometry_key': self.edges_geometry_key,
                'crs': self.crs}

    def set_nodes_coordinates(self, coordinates: dict):
        """Set nodes coordinates with a given dictionary of coordinates (can be used for a subset of all nodes).

        Parameters
        ----------
        coordinates: dict :
            Dictionary mapping node names and two-element list of coordinates.
        """
        for n, coords in coordinates.items():
            node_data = self.nodes[n]
            node_data[self.nodes_geometry_key] = Point(coords)

    def to_nx_class(self):
        """Return the closest networkx class (in the inheritance graph)."""
        return nx.Graph

    def copy(self, as_view=False):
        """Return a copy of the graph (see ``networkx.Graph.copy``)."""
        nx_graph_class = self.to_nx_class()
        graph = nx_graph_class.copy(self, as_view)
        return self.__class__(graph)

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph (see ``networkx.Graph.to_directed``)."""
        if as_view:
            return nx.Graph.to_directed(self, as_view)
        else:
            graph_class = self.to_directed_class()
            directed_graph = nx.Graph.to_directed(self, as_view)
            return graph_class(directed_graph)

    def to_directed_class(self):
        """Returns the class to use for empty directed copies (see ``networkx.Graph.to_directed_class``)."""
        return gnx.GeoDiGraph

    def to_undirected(self, as_view=False):
        """Return an undirected copy of the graph (see ``networkx.Graph.to_undirected``)."""
        if as_view:
            return nx.Graph.to_undirected(self, as_view)
        else:
            graph_class = self.to_undirected_class()
            undirected_graph = nx.Graph.to_undirected(self, as_view)
            return graph_class(undirected_graph)

    def to_undirected_class(self):
        """Returns the class to use for empty undirected copies (see ``networkx.Graph.to_undirected_class``)."""
        return gnx.GeoGraph

    def to_crs(self, crs=None, epsg=None, inplace=False):
        """Transform nodes and edges geometries to a new coordinate reference system.

        Parameters
        ----------
        crs : dict or str
             Output projection parameters as string or in dictionary form (Default value = None).
        epsg : int
             EPSG code specifying output projection.
        inplace : bool
             If True, the modification is done inplace, otherwise a new graph is created (Default value = False).

        Returns
        -------
        None or GeoGraph
            Nothing is returned if the transformation is inplace, a new GeoGraph is returned otherwise.

        See Also
        --------
        geopandas.GeoSeries.to_crs
        """
        if self.crs is None:
            raise ValueError('Cannot transform naive geometries. Please set a crs on the graph first.')
        if inplace:
            graph = self
        else:
            graph = self.copy()
        # Get new nodes coordinates
        nodes_as_points = graph.get_nodes_as_point_series()
        transformed_nodes = nodes_as_points.to_crs(crs, epsg)
        # Get new edges coordinates
        edges_as_lines = graph.get_edges_as_line_series()
        transformed_edges = edges_as_lines.to_crs(crs, epsg)
        # Operate the transformation
        for n, point in transformed_nodes.iteritems():
            node_data = graph.nodes[n]
            node_data[graph.nodes_geometry_key] = point
        for e, line in transformed_edges.iteritems():
            edge_data = graph.edges[e]
            edge_data[graph.edges_geometry_key] = line
        graph.crs = transformed_nodes.crs
        if not inplace:
            return graph

    def add_node(self, node_for_adding, geometry=None, **attr):
        """Add a single node `node_for_adding` with its given geometry.

        See Also
        --------
        nx.Graph.add_node add_nodes_from

        Examples
        --------
        >>> import geonetworkx as gnx
        >>> g = gnx.GeoGraph()
        >>> g.add_node(1, gnx.Point(2, 3))
        >>> print(g.nodes[1]["geometry"])
        POINT (2 3)
        """
        if geometry is not None:
            attr.update({self.nodes_geometry_key: geometry})
        super().add_node(node_for_adding, **attr)

    def add_nodes_from(self, nodes_for_adding, **attr):
        """Add multiple nodes with potentially given geometries.

        If no geometry is provided, behaviour is same as the
        ``nx.Graph.add_nodes_from`` method.

        See Also
        --------
        nx.Graph.add_nodes_from add_node

        Examples
        --------
        >>> import geonetworkx as gnx
        >>> g = gnx.GeoGraph()
        >>> g.add_nodes_from([(1, gnx.Point(1, 1)),
        ...                   (2, gnx.Point(2, 1)),
        ...                   (3, gnx.Point(3, 1))])
        >>> print(g.nodes[2]["geometry"])
        POINT (2 1)
        """
        geom_key = self.nodes_geometry_key

        def format_element(n):
            try:
                nn, p = n
                if isinstance(p, Point):
                    return nn, {geom_key: p}
            except (TypeError, ValueError):
                pass
            return n
        nodes_for_adding_formatted = map(format_element, nodes_for_adding)
        super().add_nodes_from(nodes_for_adding_formatted, **attr)

    def _get_nodes_geometries_from_edge_geometry(self, u, v, geometry):
        """For each node of the edge, return the node geometry deduced from the linestring if it not already present."""
        u_geometry = v_geometry = None
        if isinstance(geometry, LineString):
            if u not in self._node:
                u_geometry = Point(geometry.coords[0])
            if v not in self._node:
                v_geometry = Point(geometry.coords[-1])
        return u_geometry, v_geometry

    def add_edge(self, u_of_edge, v_of_edge, **attr):
        """Add an edge between u and v.

        If one of the node is not already in the graph and a geometry is provided, the node geometry is deduced from
        the first or last point of the linestring.

        Examples
        --------
        >>> import geonetworkx as gnx
        >>> g = gnx.GeoGraph()
        >>> g.add_edge(1, 2, geometry=gnx.LineString([(0, 0), (1, 1)]))
        >>> print(g.nodes[2]["geometry"])
        POINT (1 1)
        """
        u_geometry, v_geometry = self._get_nodes_geometries_from_edge_geometry(u_of_edge, v_of_edge,
                                                                               attr.get(self.edges_geometry_key, None))
        self.to_nx_class().add_edge(self, u_of_edge, v_of_edge, **attr)
        if u_geometry is not None:
            self.nodes[u_of_edge][self.nodes_geometry_key] = u_geometry
        if v_geometry is not None:
            self.nodes[v_of_edge][self.nodes_geometry_key] = v_geometry

    def _get_nodes_geometries_to_set_for_edges_adding(self, ebunch_to_add, attr):
        """Return a dictionary of nodes geometries to set when adding a set of edges."""
        nodes_geometry_to_set = dict()
        for e in ebunch_to_add:
            ne = len(e)
            edge_geometry = attr.get(self.edges_geometry_key, None)
            if ne == 3:
                try:
                    edge_geometry = e[2].get(self.edges_geometry_key, None)
                except:
                    pass
            elif ne == 4 and self.edges_geometry_key in e[3]:
                edge_geometry = e[3][self.edges_geometry_key]
            u, v = e[0], e[1]
            u_geometry, v_geometry = self._get_nodes_geometries_from_edge_geometry(u, v, edge_geometry)
            if u_geometry is not None and u not in nodes_geometry_to_set:
                nodes_geometry_to_set[u] = u_geometry
            if v_geometry is not None and v not in nodes_geometry_to_set:
                nodes_geometry_to_set[v] = v_geometry
        return nodes_geometry_to_set

    def add_edges_from(self, ebunch_to_add, **attr):
        """Add all the edges in ebunch_to_add and add nodes geometry if they are not present.

        If one of the node is not already in the graph and a geometry is provided, the node geometry is deduced from
        the first or last point of the linestring.

        Examples
        --------
        >>> import geonetworkx as gnx
        >>> g = gnx.GeoGraph()
        >>> g.add_edges_from([(0, 1, dict(geometry=gnx.LineString([(0, 0), (1, 1)]))),
        ...                  (1, 2, dict(geometry=gnx.LineString([(1, 1), (2, 2)])))])
        >>> print(g.nodes[2]["geometry"])
        POINT (2 2)

        >>> g = gnx.GeoMultiGraph()
        >>> g.add_edges_from([(0, 1, 7, dict(geometry=gnx.LineString([(-1, 0), (1, 1)]))),
        ...                   (1, 2, 8, dict(geometry=gnx.LineString([(1, 1), (2, 2)])))])
        [7, 8]
        >>> print(g.nodes[1]["geometry"])
        POINT (1 1)

        See Also
        --------
        add_edge
        nx.Graph.add_edges_from

        """
        r1, r2 = itertools.tee(ebunch_to_add)  # in case a generator is passed
        nodes_geometry_to_set = self._get_nodes_geometries_to_set_for_edges_adding(r1, attr)
        result = super().add_edges_from(r2, **attr)
        for u, g in nodes_geometry_to_set.items():
            self.nodes[u][self.nodes_geometry_key] = g
        if result is not None:
            return result

    def to_utm(self, inplace=False):
        """Project graph coordinates to the corresponding UTM (Universal Transverse Mercator)

        Parameters
        ----------
        inplace : bool
             If True, the modification is done inplace, otherwise a new graph is returned (Default value = False).

        Example
        -------
        >>> import geonetworkx as gnx
        >>> from shapely.geometry import Point
        >>> g = gnx.GeoGraph(crs=gnx.WGS84_CRS)
        >>> g.add_edge(1, 2, geometry=gnx.LineString([(4.28, 45.5), (4.32, 45.48)]))
        >>> g.to_utm(inplace=True)
        >>> print(g.crs)
        +proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +type=crs
        >>> print(g.nodes[1]["geometry"])
        POINT (600002.1723317318 5039293.296216004)

        See Also
        --------
        to_crs
        """
        graph_centroid = self.get_nodes_as_multipoint().centroid
        utm_crs = gnx.get_utm_crs(graph_centroid)
        if inplace:
            self.to_crs(utm_crs, inplace=True)
        else:
            return self.to_crs(utm_crs, inplace=False)

    def nodes_to_gdf(self) -> gpd.GeoDataFrame:
        """Create a ``geopandas.GeoDataFrame`` from nodes of the current graph. The column representing the geometry is
        named after the current ``nodes_geometry_key`` attribute.

        Returns
        -------
        gpd.GeoDataFrame
            The resulting GeoDataFrame : one row is a node

        See Also
        --------
        get_nodes_as_point_series, edges_to_gdf

        """
        nodes = {node: data for node, data in self.nodes(data=True)}
        gdf_nodes = gpd.GeoDataFrame(nodes).T
        gdf_nodes[settings.NODE_ID_COLUMN_NAME] = gdf_nodes.index
        gdf_nodes.set_geometry(self.nodes_geometry_key, inplace=True)
        gdf_nodes.crs = self.crs
        return gdf_nodes

    def edges_to_gdf(self) -> gpd.GeoDataFrame:
        """Create a ``gpd.GeoDataFrame`` from edges of the current graph. The column representing the geometry is
        named after the current ``edges_geometry_key`` attribute.

        Returns
        -------
        gdf_edges: geopandas.GeoDataFrame
            The resulting GeoDataFrame : one row is an edge

        See Also
        --------
        get_edges_as_line_series, nodes_to_gdf

        """
        # create a list to hold our edges, then loop through each edge in the graph
        edges = []
        for u, v, data in self.edges(data=True):
            # for each edge, add key and all attributes in data dict to the edge_details
            edge_details = {settings.EDGE_FIRST_NODE_COLUMN_NAME: u, settings.EDGE_SECOND_NODE_COLUMN_NAME: v}
            edge_details.update(data)
            # if edge doesn't already have a geometry attribute, create one now
            if self.edges_geometry_key not in data:
                point_u = self.nodes[u][self.nodes_geometry_key]
                point_v = self.nodes[v][self.nodes_geometry_key]
                edge_details[self.edges_geometry_key] = LineString([point_u, point_v])
            edges.append(edge_details)
        # create a GeoDataFrame from the list of edges and set the CRS
        gdf_edges = gpd.GeoDataFrame(edges)
        gdf_edges.set_geometry(self.edges_geometry_key, inplace=True)
        gdf_edges.crs = self.crs
        return gdf_edges

    def add_nodes_from_gdf(self, gdf: gpd.GeoDataFrame, node_index_attr=None):
        """Add nodes with the given `GeoDataFrame` and fill nodes attributes with the geodataframe columns.

        Parameters
        ----------
        gdf :
            GeoDataFrame representing nodes to add (one row for one node).
        node_index_attr :
            Node index attribute for labeling nodes. If ``None``, the dataframe index is used, else
            the given column is used. (Default value = None)

        See Also
        --------
        add_edges_from_gdf

        """
        if not (gnx.is_null_crs(self.crs) or gnx.is_null_crs(gdf.crs) or gnx.crs_equals(gdf.crs, self.crs)):
            gdf = gdf.to_crs(self.crs, inplace=False)
        if node_index_attr is not None:
            gdf = gdf.set_index(node_index_attr, drop=True, inplace=False)
        if gdf._geometry_column_name != self.nodes_geometry_key:
            gdf = gdf.rename(columns={gdf._geometry_column_name: self.nodes_geometry_key}, inplace=False)
            gdf.set_geometry(self.nodes_geometry_key, inplace=True)
        self.add_nodes_from(gdf.iterrows())
        self.check_nodes_validity()

    def add_edges_from_gdf(self, gdf: gpd.GeoDataFrame, edge_first_node_attr=None, edge_second_node_attr=None):
        """Add edges with the given `GeoDataFrame`. If no dataframe columns are specified for first and second node,
        the dataframe index must be a multi-index `(u, v)`.

        Parameters
        ----------
        gdf :
            GeoDataFrame representing edges to add (one row for one edge).
        edge_first_node_attr :
            Edge first node attribute. If ``None``, the dataframe index is used, else the given
            column is used. Must be used with ``edge_second_node_attr``. (Default value = None)
        edge_second_node_attr :
            Edge second node attribute. If ``None``, the dataframe index is used, else the
            given column is used. Must be used with ``edge_first_node_attr``. (Default value = None)

        See Also
        --------
        add_nodes_from_gdf
        """
        if not (gnx.is_null_crs(self.crs) or gnx.is_null_crs(gdf.crs) or gnx.crs_equals(gdf.crs, self.crs)):
            gdf = gdf.to_crs(self.crs, inplace=False)
        if edge_first_node_attr is not None and edge_second_node_attr is not None:
            gdf = gdf.set_index([edge_first_node_attr, edge_second_node_attr], drop=True, inplace=False)
        if gdf._geometry_column_name != self.edges_geometry_key:
            gdf = gdf.rename(columns={gdf._geometry_column_name: self.edges_geometry_key}, inplace=False)
            gdf.set_geometry(self.edges_geometry_key, inplace=True)
        self.add_edges_from((*r[0], r[1]) for r in gdf.iterrows())
