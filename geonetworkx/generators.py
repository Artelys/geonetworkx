# -*- coding: utf-8 -*-
import networkx as nx
from geonetworkx.geograph import GeoGraph
from geonetworkx import settings
from geonetworkx.geometry_operations import split_line
from geonetworkx.utils.geograph_utils import get_line_start


def _get_ego_boundaries(graph: GeoGraph, ego_graph: GeoGraph, n, radius: float, distance=None) -> tuple:
    """Retrieve all information to build an extended ego-graph. See ``gnx.extended_ego_graph`` for more info."""
    if graph.is_multigraph():
        out_edges_options = {"keys": True}
    else:
        out_edges_options = {}
    boundary_nodes = []
    boundary_edges = dict()
    for u in ego_graph.nodes:
        node_edges = graph.edges(u, **out_edges_options)
        for e in node_edges:
            v = e[1]
            if v in ego_graph.nodes:
                continue
            edge_data = graph.edges[e]
            if graph.edges_geometry_key not in edge_data:
                continue
            # The distance d(n, u) has already been computed in the networkx ego graph method,
            # but for the sake of code simplicity, we are doing it again here.
            d_n_u = nx.dijkstra_path_length(graph, n, u, weight=distance)
            edge_distance = edge_data.get(distance, 1)
            distance_margin = radius - d_n_u
            interpolation_factor = distance_margin / edge_distance
            edge_geometry = edge_data[graph.edges_geometry_key]
            edge_length = edge_geometry.length
            if get_line_start(graph, e, edge_geometry) == u:
                interpolation_distance = interpolation_factor * edge_length
                inside_line_part_index = 0
            else:
                interpolation_distance = (1.0 - interpolation_factor) * edge_length
                inside_line_part_index = 1
            if interpolation_distance <= 0 or interpolation_distance >= edge_length:
                continue
            b_node_point = edge_geometry.interpolate(interpolation_distance)
            b_node_name = settings.BOUNDARY_NODE_PREFIX + str(u) + '_' + str(v)
            boundary_nodes.append((b_node_name, {graph.nodes_geometry_key: b_node_point}))
            split_edge = split_line(edge_geometry, interpolation_distance)
            inside_edge_data = {graph.edges_geometry_key: split_edge[inside_line_part_index],
                                distance: distance_margin}
            outside_edge_data = {graph.edges_geometry_key: split_edge[1 - inside_line_part_index],
                                 distance: edge_distance - distance_margin}
            boundary_edges[e] = [(u, b_node_name, inside_edge_data), (b_node_name, v, outside_edge_data)]
    return boundary_nodes, boundary_edges


def extended_ego_graph(graph: GeoGraph, n, radius=1, center=True, undirected=False, distance=None) -> GeoGraph:
    """Returns induced subgraph of neighbors centered at node n within a given radius extended by interpolated nodes on
    boundary edges.

    Note that the returned graph is not a subgraph of the input graph because it will have boundary nodes in addition.
    It means that a node is added on each edge leaving the ego-graph to represent the furthest reachable point on
    this edge. The boundary node is added at given node using a linear interpolation. A boundary node :math:`b` will
    be added on the edge :math:`(u, v)` if :math:`d(n, u) < radius < d(n, v)`. The boundary will be placed along the
    edge geometry at the following distance:

    .. math::
        d(u, b) =  \\frac{\\text{radius} - d(n, u)}{d(u, v)}

    Furthermore, the attribute ``distance`` is filled with the value :math:`d(u, b)`.

    Parameters
    ----------
    graph : GeoGraph, GeoDiGraph, GeoMultiGraph or GeoMultiDiGraph
        A Geograph or subclass
    n :
        A single source node
    radius : float or int
        Include all neighbors of distance<=radius from n. (Default value = 1)
    center : bool
        If False, do not include center node in graph (Default value = True)
    undirected : bool
        If True use both in- and out-neighbors of directed graphs. (Default value = False)
    distance : str
        Use specified edge data key as distance.  For example, setting distance='weight' will use the edge
        weight to measure the distance from the node n. (Default value = None)

    Returns
    -------
    GeoGraph, GeoDiGraph, GeoMultiGraph or GeoMultiDiGraph
        The resulting graph with node, edge, and graph attributes copied.

    See Also
    --------
    add_ego_boundary_nodes

    """
    ego_graph = nx.ego_graph(graph, n, radius, center=True, undirected=undirected, distance=distance)
    if undirected:
        working_graph = graph.to_undirected()
    else:
        working_graph = graph
    # Retrieve boundary edge and cut geometries
    boundary_nodes, boundary_edges = _get_ego_boundaries(working_graph, ego_graph, n, radius, distance)
    ego_graph.add_nodes_from(boundary_nodes)
    edges_to_add = []
    for e, ego_edges in boundary_edges.items():
        edges_to_add.append(ego_edges[0])
    ego_graph.add_edges_from(edges_to_add)
    if not center:
        ego_graph.remove_node(n)
    return ego_graph


def add_ego_boundary_nodes(graph: GeoGraph, n, radius, distance, undirected=False):
    """Modify the given graph to add boundary nodes at exact radius distance. An edge :math:`(u, v)` is a boundary edge
    if :math:`(u, v)` if :math:`d(n, u) < radius < d(n, v)`. A boundary node is added on the edge to represent the ego-
    graph limit. See ``gnx.extended_ego_graph`` for more info.

    Parameters
    ----------
    graph : GeoGraph, GeoDiGraph, GeoMultiGraph or GeoMultiDiGraph
        Input graph to modify
    n :
        A single source node
    radius : float or int
        Include all neighbors of distance<=radius from n.
    distance : str
        Use specified edge data key as distance. For example, setting distance='weight' will use the edge
        weight to measure the distance from the node n.
    undirected : bool
        If True use both in- and out-neighbors of directed graphs. (Default value = False)

    See Also
    --------
    extended_ego_graph

    """
    ego_graph = nx.ego_graph(graph, n, radius, center=True, undirected=undirected, distance=distance)
    if undirected:
        working_graph = graph.to_undirected()
    else:
        working_graph = graph
    # Retrieve boundary edge and cut geometries
    boundary_nodes, boundary_edges = _get_ego_boundaries(working_graph, ego_graph, n, radius, distance)
    graph.add_nodes_from(boundary_nodes)
    edges_to_add = []
    for e, ego_edges in boundary_edges.items():
        edges_to_add.extend(ego_edges)
        graph.remove_edge(*e)
    graph.add_edges_from(edges_to_add)
