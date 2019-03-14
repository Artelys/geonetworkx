"""
    File name: tools
    Author: Artelys
    Creation date: 08/01/2019
    Python Version: 3.6
"""
import numpy as np
import networkx as nx
from networkx.classes.filters import no_filter
import geopandas as gpd
from shapely.geometry import Point, LineString
from geonetworkx.geograph import GeoGraph
import geonetworkx.settings as settings
from geonetworkx.geometry_operations import get_closest_line_from_points, split_line, coordinates_almost_equal
from geonetworkx.utils import get_new_node_unique_name, euclidian_distance, get_line_ordered_edge, is_nan
from collections import defaultdict


def spatial_points_merge(graph: GeoGraph, points_gdf: gpd.GeoDataFrame, inplace=False, merge_direction="both",
                         node_filter=no_filter, edge_filter=no_filter, intersection_nodes_attr=None) -> GeoGraph:
    """
    Merge given points as node with a spatial merge. Points are projected on the closest edge of the
    graph and an intersection node is added if necessary. If two nodes a given point and a node have the same name, with
    equal coordinates, then the node is considered as already in the graph. A discretization tolerance
    (``settings.DISCRETIZATION_TOLERANCE``) is used for edges lines and is set by default to a constant matching the
    WGS84 crs. If another crs is used, results may be inconsistent (high computational time or inaccuracy). New nodes
    created from the geodataframe have attributes described by other columns (except if an attribute value is `nan`).

    :param graph: A GeoGraph or derived class describing a spatial graph.
    :param points_gdf: A list of point describing new nodes to add.
    :param inplace: If True, do operation inplace and return None.
    :param merge_direction: For directed graphs only:

         * ``'both'``: 2 edges are added: graph -> new node and new node -> graph
         * ``'in'``: 1 edge is added: new_node -> graph
         * ``'out'``: 1 edge is added: graph -> new_node
    :param node_filter: A node filter (lambda) to exclude nodes (and by the way all concerned edges) from the projection
     operation.
    :param edge_filter: An edge filter (lambda) to exclude edges on which the projection will not take place.
    :param intersection_nodes_attr: A dictionary of attributes (constant for all added intersection nodes).
    :return: None if inplace, new graph otherwise.
    """
    if not inplace:
        graph = graph.copy()
    # 1. Find closest edge for each point
    graph_view = nx.graphviews.subgraph_view(graph, filter_node=node_filter, filter_edge=edge_filter)
    edges_as_lines = nx.get_edge_attributes(graph_view, graph.edges_geometry_key)
    if len(edges_as_lines) == 0:
        raise ValueError("No edge geometry has been found in the given merging edges, at least one edge geometry is"
                         " required for a merge operation")
    points = points_gdf.geometry
    points_coords = np.array([[p.x, p.y] for p in points])
    lines_indexes = get_closest_line_from_points(points_coords, edges_as_lines.values())
    edges_to_split = defaultdict(dict)
    # Add node, intersection node and edge (node, intersection node)
    for p, p_index, point in zip(range(len(points_gdf)), points_gdf.index, points):
        # 1.1 Add given node
        if p_index in graph.nodes:
            if coordinates_almost_equal([point.x, point.y], graph.get_node_coordinates(p_index)):
                continue
            else:
                node_name = get_new_node_unique_name(graph, p_index)
        else:
            node_name = p_index
        node_info = {c: points_gdf.at[p_index, c] for c in points_gdf.columns if not is_nan(points_gdf.at[p_index, c])}
        node_info[graph.nodes_geometry_key] = point
        graph.add_node(node_name, **node_info)
        # 1.2 Add projected node if necessary
        closest_edge_name = list(edges_as_lines.keys())[lines_indexes[p]]
        closest_line = edges_as_lines[closest_edge_name]
        closest_line_length = closest_line.length
        intersection_distance_on_line = closest_line.project(point)
        # if the intersection point is on the edge
        if 0 < intersection_distance_on_line < closest_line_length:
            projected_point = closest_line.interpolate(intersection_distance_on_line)
            intersection_node_name = get_new_node_unique_name(graph, settings.INTERSECTION_PREFIX + str(p_index))
            intersection_node_info = {graph.nodes_geometry_key: projected_point}
            if intersection_nodes_attr is not None:
                intersection_node_info.update(intersection_nodes_attr)
            graph.add_node(intersection_node_name, **intersection_node_info)
            # Store line to modify
            edges_to_split[closest_edge_name][intersection_node_name] = intersection_distance_on_line
        else:  # if the intersection point is on of the two edge extremities
            first_node = closest_edge_name[0]
            first_node_point = graph.nodes[first_node][graph.nodes_geometry_key]
            second_node = closest_edge_name[1]
            second_node_point = graph.nodes[second_node][graph.nodes_geometry_key]
            distance_to_first_extremity = euclidian_distance(point, first_node_point)
            distance_to_second_extremity = euclidian_distance(point, second_node_point)
            if distance_to_first_extremity < distance_to_second_extremity:
                intersection_node_name = first_node
            else:
                intersection_node_name = second_node
        # 1.3 Add edge : node <-> intersection_node
        in_edge_data = {graph.edges_geometry_key: LineString([graph.get_node_coordinates(node_name),
                                                              graph.get_node_coordinates(intersection_node_name)])}
        if graph.is_directed():
            out_edge_data = {graph.edges_geometry_key: LineString([graph.get_node_coordinates(intersection_node_name),
                                                                   graph.get_node_coordinates(node_name)])}
            if merge_direction == "both":
                graph.add_edge(node_name, intersection_node_name, **in_edge_data)
                graph.add_edge(intersection_node_name, node_name, **out_edge_data)
            elif merge_direction == "in":
                graph.add_edge(node_name, intersection_node_name, **in_edge_data)
            else:  # "out"
                graph.add_edge(intersection_node_name, node_name, **out_edge_data)
        else:
            graph.add_edge(node_name, intersection_node_name, **in_edge_data)
    # 2. Split edges where a node have been projected
    for e in edges_to_split:
        intersection_nodes = edges_to_split[e]
        if len(intersection_nodes) > 0:
            initial_line = edges_as_lines[e]
            # 2.1 remove initial edge
            if graph.has_edge(*e):
                graph.remove_edge(*e)
            # 2.2 cut the initial line
            sorted_intersection_nodes = sorted(intersection_nodes.keys(), key=lambda n: intersection_nodes[n])
            distances_on_initial_line = [intersection_nodes[n] for n in sorted_intersection_nodes]
            split_lines = []
            cut_lines = split_line(initial_line, distances_on_initial_line[0])
            split_lines.append(cut_lines[0])
            for i in range(len(sorted_intersection_nodes) - 1):
                cut_lines = split_line(cut_lines[1], distances_on_initial_line[i + 1] - distances_on_initial_line[i])
                split_lines.append(cut_lines[0])
            split_lines.append(cut_lines[1])
            # 2.2 add intermediary edges
            oriented_edge = get_line_ordered_edge(graph, e, initial_line)
            first_edge_data = {graph.edges_geometry_key: split_lines[0]}
            graph.add_edge(oriented_edge[0], sorted_intersection_nodes[0], **first_edge_data)
            last_edge_data = {graph.edges_geometry_key: split_lines[-1]}
            graph.add_edge(sorted_intersection_nodes[-1], oriented_edge[1], **last_edge_data)
            for i in range(len(sorted_intersection_nodes) - 1):
                edge_data = {graph.edges_geometry_key: split_lines[i + 1]}
                graph.add_edge(sorted_intersection_nodes[i], sorted_intersection_nodes[i + 1], **edge_data)
    if not inplace:
        return graph


def spatial_graph_merge(base_graph: GeoGraph, other_graph: GeoGraph,
                        inplace=False, merge_direction="both", node_filter=None, intersection_nodes_attr=None):
    """
    Operates spatial merge between two graphs. Spatial edge projection is used on merging nodes (see
    ``spatial_points_merge``). The ``base_graph`` attributes have higher priority than the ``other_graph`` attributes (
    i.e. if graphs have common graph attributes, nodes or edges, the ``base_graph`` attributes will be kept).

    :param base_graph: Base graph on which the merge operation is done.
    :param other_graph: Input graph to merge. Modified graph if operation is done inplace.
    :param inplace: If True, do operation inplace and return None.
    :param merge_direction: See ``spatial_points_merge``
    :param node_filter: Lambda returning if a given node (from the ``other_graph`` graph) has to be merged.
    :param intersection_nodes_attr: A dictionary of attributes (constant for all added intersection nodes).
    :return: A new graph with the same type as ``base_graph`` if not inplace.
    """
    if base_graph.is_directed() != other_graph.is_directed():
        raise ValueError("Merging a directed graph and an undirected graph is ambiguous")
    if base_graph.is_multigraph() != other_graph.is_multigraph():
        raise ValueError("Given graphs must both be graphs or multigraph")
    if node_filter is not None:
        other_graph_view = nx.graphviews.subgraph_view(other_graph, node_filter)
    else:
        other_graph_view = other_graph
    if other_graph_view.number_of_nodes() == 0:
        raise ValueError("Given graph has no nodes to project.")
    nodes_gdf = other_graph_view.nodes_to_gdf()
    if inplace:
        spatial_points_merge(base_graph, nodes_gdf, inplace=inplace, merge_direction=merge_direction,
                             intersection_nodes_attr=intersection_nodes_attr)
        merged_graph = base_graph
    else:
        merged_graph = spatial_points_merge(base_graph, nodes_gdf, inplace=inplace, merge_direction=merge_direction,
                                            intersection_nodes_attr=intersection_nodes_attr)
    merged_graph = nx.compose(other_graph, merged_graph)
    if not inplace:
        return merged_graph


