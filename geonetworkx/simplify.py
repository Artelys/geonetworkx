import numpy as np
import networkx as nx
from shapely.geometry import Polygon, MultiPolygon
from geonetworkx import GeoGraph
from geonetworkx.utils import is_nan
from typing import Union


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

    def trim_data(data):
        keys_to_remove = set()
        for k, v in data.items():
            if remove_none and v is None:
                keys_to_remove.add(k)
            if remove_nan and is_nan(v):
                keys_to_remove.add(k)
        for k in keys_to_remove:
            del data[k]

    for n, data in used_graph.nodes(data=True):
        trim_data(data)
    for u, v, data in used_graph.edges(data=True):
        trim_data(data)
    if copy:
        return used_graph


