# -*- coding: utf-8 -*-
import networkx as nx
from shapely.geometry import Polygon, MultiPolygon, LineString
from geonetworkx import GeoGraph, GeoDiGraph
from geonetworkx.utils import is_nan
from geonetworkx.geometry_operations import merge_two_lines_with_closest_extremities
from typing import Union
from networkx.classes.filters import no_filter


def remove_isolates(graph: nx.Graph) -> int:
    """
    Removes all isolates nodes in the given graph.

    :param graph: A graph on which to remove all isolates
    :return: Number of removed isolates
    """
    isolates = list(nx.isolates(graph))
    graph.remove_nodes_from(isolates)
    return len(isolates)


def remove_self_loop_edges(graph: nx.Graph) -> int:
    """
    Remove self loop edges on nodes of the given graph.

    :param graph: A graph on which to remove all self loops.
    :return: The number of removed self loops
    """
    self_loops_edges = list(nx.selfloop_edges(graph))
    graph.remove_edges_from(self_loops_edges)
    return len(self_loops_edges)


def remove_small_connected_components(graph: nx.Graph, minimum_allowed_size: int) -> int:
    """
    Remove all connected components having strictly less than ``minimum_allowed_size``.

    :param graph: The graph on which to remove connected components
    :param minimum_allowed_size: The minimum number of nodes where a connected component is kept.
    :return: The number of removed connected components
    """
    connected_components = list(nx.connected_components(graph))
    nb_removed_cc = 0
    for c_ix, cc in enumerate(connected_components):
        if len(cc) < minimum_allowed_size:
            graph.remove_nodes_from(cc)
            nb_removed_cc += 1
        else:
            for n in cc:
                graph.nodes[n]['cc'] = c_ix
    return nb_removed_cc


def trim_graph_with_polygon(graph: GeoGraph, polygon: Union[Polygon, MultiPolygon], copy=False, method="intersects"):
    """
    Trim a graph with a given polygon. Keep only the nodes and edges that intersect (or are within) the polygon.

    :param graph: A GeoGraph (or subclass)
    :param polygon: A ``shapely.Polygon`` describing the area to keep
    :param copy: If ``True``, a deep copy is done and a new graph is returned.
    :param method: If set to ``"intersects"``, the ``shapely.intersects`` is used (keeps nodes and edges that
        intersects the polygon). If set to ``"within"``, the ``shapely.within`` is used (keep nodes and edges that are
        strictly into the polygon).
    :return: The modified graph if ``copy`` is ``True``.
    """
    if copy:
        used_graph = graph.copy()
    else:
        used_graph = graph
    if method not in ["intersects", "within"]:
        raise ValueError("Unknown method for trimming : '%s'" % str(method))
    nodes_series = used_graph.get_nodes_as_point_series()
    edges_as_series = used_graph.get_edges_as_line_series()
    if method == 'intersects':
        nodes_criteria = ~ nodes_series.intersects(polygon)
        edges_criteria = ~ edges_as_series.intersects(polygon)
    else:
        nodes_criteria = ~ nodes_series.within(polygon)
        edges_criteria = ~ edges_as_series.within(polygon)
    nodes_to_remove = nodes_series[nodes_criteria].index
    edges_to_remove = edges_as_series[edges_criteria].index
    used_graph.remove_nodes_from(nodes_to_remove)
    used_graph.remove_edges_from(edges_to_remove)
    if copy:
        return used_graph


def remove_nan_attributes(graph: nx.Graph, remove_nan=True, remove_none=True, copy=False):
    """
    Remove the `nan` and `None` values from nodes and edges attributes.

    :param graph: Graph (or subclass)
    :param remove_nan: If true, remove the `nan` values (test is ``val is np.nan``)
    :param remove_none: If true, remove the ``None`` values (test is ``val is None``)
    :param copy: If True, a copy of the graph is returned, otherwise the graph is modified inplace.
    :return: The modified graph if ``copy`` is true.
    """
    if copy:
        used_graph = graph.copy()
    else:
        used_graph = graph

    def trim_data(d):
        keys_to_remove = set()
        for k, val in d.items():
            if remove_none and val is None:
                keys_to_remove.add(k)
            if remove_nan and is_nan(val):
                keys_to_remove.add(k)
        for k in keys_to_remove:
            del d[k]

    for n, data in used_graph.nodes(data=True):
        trim_data(data)
    for u, v, data in used_graph.edges(data=True):
        trim_data(data)
    if copy:
        return used_graph


def get_dead_ends(graph: nx.Graph, node_filter=no_filter, only_strict=False):
    """Return the list of dead end in the given graph. A dead end is defined as a node having only one neighbor. For
    directed graphs, a strict dead end is a node having a unique predecessor and no successors. A weak dead end is a
    node having a unique predecessor that is also its unique successor.

    :param graph: Graph to parse.
    :param node_filter: Evaluates to true if a node can be considered as dead end, false otherwise.
    :param only_strict: If true, remove only strict dead ends. Used only for directed graphs.
    """
    if graph.is_directed():
        dead_ends = []
        for n in graph.nodes():
            if not node_filter(n):
                continue
            nb_predecessors = len(graph.pred[n])
            if nb_predecessors != 1:
                continue
            nb_successors = len(graph.succ[n])
            if nb_successors == 0:
                dead_ends.append(n)
            elif not only_strict and nb_successors == 1:
                pred = next(iter(graph.pred[n]))
                if pred in graph.successors(n):
                    dead_ends.append(n)
        return dead_ends
    else:
        return [n for n in graph.nodes() if node_filter(n) and len(graph.neighbors(n)) == 1]


def remove_dead_ends(graph: nx.Graph, node_filter=no_filter, only_strict=False):
    """Remove dead ends from a given graph. A dead end is defined as a node having only one neighbor. For
    directed graphs, a strict dead end is a node having a unique predecessor and no successors. A weak dead end is a
    node having a unique predecessor that is also its unique successor.

    :param graph: Graph to simplify
    :param node_filter: Evaluates to true if a node can be removed, false otherwise.
    :param only_strict: If true, remove only strict dead ends. Used only for directed graphs.
    """
    nodes_to_remove = get_dead_ends(graph, node_filter, only_strict)
    while nodes_to_remove:
        graph.remove_nodes_from(nodes_to_remove)
        nodes_to_remove = get_dead_ends(graph, node_filter, only_strict)


def _clean_merge_mapping(edge_mapping: dict, new_edge: tuple, old_edges: list, directed: bool):
    """For the two-degree node merge operation, it cleans the new-old edges mapping dictionary by reporting original
    edges to the newest edge. It makes sure that all edges in the mapping dictionary dict are in the resulting graph."""
    for e in old_edges:
        old_edge = None
        if e in edge_mapping.keys():
            old_edge = e
        elif not directed:
            reversed_edge = (e[1], e[0], *e[2:])
            if reversed_edge in edge_mapping.keys():
                old_edge = reversed_edge
        if old_edge is not None:
            if e in edge_mapping[new_edge]:
                edge_mapping[new_edge].remove(e)
            edge_mapping[new_edge].extend(edge_mapping[old_edge])
            del edge_mapping[old_edge]


def two_degree_node_merge_for_directed_graphs(graph: GeoDiGraph, node_filter=no_filter) -> dict:
    """Merge edges that connects two nodes with a unique third node. A potential node to merge `n` must have exactly two
    different neighbors `u` and `v` with one of the following set of edges:
        * `(u, n)` and `(n, v)`
        * `(u, n)`, `(n, u)`, `(n, v)` and `(v, n)`

    For the first case, a merging edge `(u, v)` is added. Under the latter, two edges `(u, v)` and `(v, u)` are added.
    The added edges will have a geometry corresponding to concatenation of the two replaced edges. If a replaced edge
    doesn't have a geometry, the added edge will not have a geometry as well. Edges geometries must be well ordered
    (first node must match with line's first extremity), otherwise lines concatenation may not be consistent (see
    ``order_well_lines``).

    :param graph: Given graph to modify
    :param node_filter: Evaluates to true if a given node can be merged.
    :return: Dictionary indicating for each new edge the merged ones.
    """
    def _get_merging_line(g: GeoDiGraph, e1: tuple, e2: tuple) -> Union[LineString, None]:
        first_edge_geometry = g.edges[e1].get(g.edges_geometry_key, None)
        second_edge_geometry = g.edges[e2].get(g.edges_geometry_key, None)
        if first_edge_geometry is not None and second_edge_geometry is not None:
            return LineString(list(first_edge_geometry.coords) + list(second_edge_geometry.coords))
        else:
            return None

    merged_edges = dict()
    nodes = list(graph.nodes)
    for n in nodes:
        if not node_filter(n):
            continue
        in_degree = graph.in_degree(n)
        out_degree = graph.out_degree(n)
        merging_edges = []
        if in_degree == out_degree == 1:
            predecessor = next(iter(graph.pred[n]))
            successor = next(iter(graph.succ[n]))
            if predecessor == successor:
                continue
            if graph.is_multigraph():
                new_edge = (predecessor, successor, graph.new_edge_key(predecessor, successor))
                edges = [(predecessor, n, next(iter(graph.adj[predecessor][n]))),
                         (n, successor, next(iter(graph.adj[n][successor])))]
            else:
                new_edge = (predecessor, successor)
                edges = [(predecessor, n), (n, successor)]
            merged_line = _get_merging_line(graph, edges[0], edges[1])
            merging_edges = [(new_edge, merged_line)]
            merged_edges[new_edge] = edges.copy()  # copy for non destructive delete in clean merger function
            _clean_merge_mapping(merged_edges, new_edge, edges, True)
        if in_degree == out_degree == 2:
            successors = list(graph.succ[n])
            if all(p in successors for p in graph.pred[n]):
                if successors[1] == successors[0]:
                    continue
                if graph.is_multigraph():
                    back_new_edge = (successors[1], successors[0], graph.new_edge_key(successors[1], successors[0]))
                    forth_new_edge = (successors[0], successors[1], graph.new_edge_key(successors[0], successors[1]))
                    back_edges = [(successors[1], n, next(iter(graph.adj[successors[1]][n]))),
                                  (n, successors[0], next(iter(graph.adj[n][successors[0]])))]
                    forth_edges = [(successors[0], n, next(iter(graph.adj[successors[0]][n]))),
                                   (n, successors[1], next(iter(graph.adj[n][successors[1]])))]
                else:
                    back_new_edge = (successors[1], successors[0])
                    forth_new_edge = (successors[0], successors[1])
                    back_edges = [(successors[1], n), (n, successors[0])]
                    forth_edges = [(successors[0], n), (n, successors[1])]
                back_merged_line = _get_merging_line(graph, back_edges[0], back_edges[1])
                forth_merged_line = _get_merging_line(graph, forth_edges[0], forth_edges[1])
                merging_edges = [(back_new_edge, back_merged_line),
                                 (forth_new_edge, forth_merged_line)]
                merged_edges[back_new_edge] = back_edges.copy()
                _clean_merge_mapping(merged_edges, back_new_edge, back_edges, True)
                merged_edges[forth_new_edge] = forth_edges.copy()
                _clean_merge_mapping(merged_edges, forth_new_edge, forth_edges, True)
        if merging_edges:
            # Remove node (and thus edges)
            graph.remove_node(n)
            # Add merging edges
            for new_edge, line in merging_edges:
                merging_edge_attributes = {}
                if line is not None:
                    merging_edge_attributes[graph.edges_geometry_key] = line
                graph.add_edge(*new_edge, **merging_edge_attributes)
    return merged_edges


def two_degree_node_merge_for_undirected_graphs(graph: GeoGraph, node_filter=no_filter) -> dict:
    """Merge edges that connects two nodes with a unique third node for undirected graphs. Potential nodes to merge are
    nodes with two edges connecting two different nodes. If a replaced edge doesn't have a geometry, the added edge will
    not have a geometry as well.

    :param graph: Graph to modify
    :param node_filter: Evaluates to true if a given node can be merged.
    :return: Dictionary indicating for each new edge the merged ones.
    """
    if graph.is_multigraph():
        keys_args = {"keys": True}
    else:
        keys_args = {}
    merged_edges = dict()
    two_degree_nodes = [n for n in graph.nodes() if graph.degree(n) == 2 and node_filter(n)]
    for n in two_degree_nodes:
        edges = list(graph.edges(n, **keys_args))
        assert (len(edges) == 2)
        first_node = edges[0][1] if edges[0][1] != n else edges[0][0]
        second_node = edges[1][1] if edges[1][1] != n else edges[1][0]
        if graph.is_multigraph():
            new_edge = (first_node, second_node, graph.new_edge_key(first_node, second_node))
        else:
            new_edge = (first_node, second_node)
        if first_node == second_node:
            continue
        first_edge_geometry = graph.edges[edges[0]].get(graph.edges_geometry_key, None)
        second_edge_geometry = graph.edges[edges[1]].get(graph.edges_geometry_key, None)
        if first_edge_geometry is not None and second_edge_geometry is not None:
            merged_line = merge_two_lines_with_closest_extremities(first_edge_geometry, second_edge_geometry)
        else:
            merged_line = None
        merged_edges[new_edge] = edges.copy()
        _clean_merge_mapping(merged_edges, new_edge, edges, False)
        # Remove node (and thus edges)
        graph.remove_node(n)
        # Add merging edge
        merging_edge_attributes = {}
        if merged_line is not None:
            merging_edge_attributes[graph.edges_geometry_key] = merged_line
        graph.add_edge(*new_edge, **merging_edge_attributes)
    return merged_edges


def two_degree_node_merge(graph: GeoGraph, node_filter=no_filter) -> dict:
    """Merge edges that connects two nodes with a unique third node. See ``two_degree_node_merge_for_directed_graphs``
    or ``two_degree_node_merge_for_undirected_graphs`` for more details.

    :param graph: Graph to modify
    :param node_filter: Evaluates to true if a given node can be merged.
    :return: Dictionary indicating for each new edge the merged ones.
    """
    if graph.is_directed():
        return two_degree_node_merge_for_directed_graphs(graph, node_filter)
    else:
        return two_degree_node_merge_for_undirected_graphs(graph, node_filter)
