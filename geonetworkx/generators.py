# -*- coding: utf-8 -*-
import networkx as nx
from geonetworkx.geograph import GeoGraph
from geonetworkx import settings
from geonetworkx.geometry_operations import split_line
from geonetworkx.utils.geograph_utils import get_line_start


def extended_ego_graph(G: GeoGraph, n, radius=1, center=True, undirected=False, distance=None) -> GeoGraph:
    """Returns induced subgraph of neighbors centered at node n within a given radius extended by interpolated nodes on
    boundary edges.

    :param G: A Geograph or subclass
    :param n: A single source node
    :param radius: Include all neighbors of distance<=radius from n.
    :param center: If False, do not include center node in graph
    :param undirected: If True use both in- and out-neighbors of directed graphs.
    :param distance: Use specified edge data key as distance.  For example, setting distance='weight' will use the edge
        weight to measure the distance from the node n.
    :return: The resulting graph with node, edge, and graph attributes copied. Note that the returned graph is not a
        subgraph of the input graph because it will have boundary nodes in addition.

    It means that a node is added on each edge leaving the ego-graph to represent the furthest
    reachable point on this edge. The boundary node is added at given node using a linear interpolation. A boundary
    node :math:`b` will be added on the edge :math:`(u, v)` if :math:`d(n, u) < radius < d(n, v)`. The boundary will be
    placed along the edge geometry at the following distance:

    .. math::
        d(u, b) =  \\frac{\\text{radius} - d(n, u)}{d(u, v)}

    Furthermore, the original edge (the one that have been interpolated) id can be retrieved in the ego graph through
    the ``settings.ORIGINAL_EDGE_KEY`` attribute that is set on boundary edges.
    """
    ego_graph = nx.ego_graph(G, n, radius, center, undirected, distance)
    if undirected:
        working_graph = G.to_undirected()
    else:
        working_graph = G
    if working_graph.is_multigraph():
        out_edges_options = {"keys": True}
    else:
        out_edges_options = {}
    nodes_to_add = []
    edges_to_add = []
    for u in ego_graph.nodes:
        node_edges = working_graph.edges(u, **out_edges_options)
        for e in node_edges:
            if e[1] in ego_graph.nodes:
                continue
            edge_data = working_graph.edges[e]
            if G.edges_geometry_key not in edge_data:
                continue
            # The distance d(n, u) has already been computed in the networkx ego graph method,
            # but for the sake of code simplicity, we are doing it again here.
            d_n_u = nx.dijkstra_path_length(working_graph, n, u, weight=distance)
            edge_cost = edge_data.get(distance, 1)
            interpolation_factor = (radius - d_n_u) / edge_cost
            edge_geometry = edge_data[G.edges_geometry_key]
            edge_length = edge_geometry.length
            if get_line_start(G, e, edge_geometry) == u:
                interpolation_distance = interpolation_factor * edge_length
                line_part_index = 0
            else:
                interpolation_distance = (1.0 - interpolation_factor) * edge_length
                line_part_index = 1
            if interpolation_distance <= 0 or interpolation_distance >= edge_length:
                continue
            b_node_point = edge_geometry.interpolate(interpolation_distance)
            b_node_name = settings.BOUNDARY_NODE_PREFIX + str(u) + '_' + str(e[1])
            nodes_to_add.append((b_node_name, {G.nodes_geometry_key: b_node_point}))
            new_edge_geometry = split_line(edge_geometry, interpolation_distance)[line_part_index]
            edges_to_add.append((u, b_node_name, {G.edges_geometry_key: new_edge_geometry,
                                                  settings.ORIGINAL_EDGE_KEY: e}))
    # Adding nodes and edges at the end
    ego_graph.add_nodes_from(nodes_to_add)
    ego_graph.add_edges_from(edges_to_add)
    return ego_graph








