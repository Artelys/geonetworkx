# -*- coding: utf-8 -*-
import math
import numpy as np
import networkx as nx
from shapely.geometry import Point, LineString, MultiPoint
import geopy.distance
import pyproj
from geonetworkx.geometry_operations import coordinates_almost_equal, insert_point_in_line
from geonetworkx.geograph import GeoGraph
import geonetworkx.settings as settings
from typing import Iterable
try:
    import srtm
except ImportError:
    srtm = None


def get_crs_as_str(crs) -> str:
    """Return the given CRS as string ``pyproj.Proj`` methods."""
    proj = pyproj.Proj(crs)
    return proj.definition_string()


def is_null_crs(crs) -> bool:
    """Test for null crs values."""
    if crs is None:
        return True
    if crs == dict():
        return True
    if crs == '':
        return True
    return False


def crs_equals(crs1, crs2) -> bool:
    """Compare CRS using ``pyproj.Proj`` objects."""
    if is_null_crs(crs1) or is_null_crs(crs1):
        return False
    return get_crs_as_str(crs1) == get_crs_as_str(crs2)


def compute_vincenty(p1, p2) -> float:
    """Returns the vincenty distance in meters given points with the format (longitude, latitude) in the WGS84
    crs."""
    return geopy.distance.distance((p1[1], p1[0]), (p2[1], p2[0])).meters


def compute_vincenty_from_points(p1: Point, p2: Point) -> float:
    return compute_vincenty([p1.x, p1.y], [p2.x, p2.y])


def approx_map_unit_factor(points_coordinates, tolerance=1e-7) -> float:
    """Compute a linear approximation of the map unit factor for 1 meter. Works only for the WGS84 CRS."""
    centroid = np.array(points_coordinates)
    lower_bound = centroid
    initial_gap = tolerance
    while compute_vincenty(centroid, centroid + initial_gap) < 1.0:
        initial_gap *= 2
    upper_bound = centroid + initial_gap
    unit_point = (lower_bound + upper_bound) / 2
    distance = compute_vincenty(centroid, unit_point)
    iteration = 0
    max_iterations = 10000
    while np.abs(distance - 1.0) > tolerance:
        if distance > 1.0:
            upper_bound = unit_point
        else:
            lower_bound = unit_point
        unit_point = (lower_bound + upper_bound) / 2
        distance = compute_vincenty(centroid, unit_point)
        iteration += 1
        if iteration > max_iterations:
            break
    return np.linalg.norm(centroid - unit_point)


def measure_line_distance(line: LineString) -> float:
    """
    Measure the length of a shapely LineString object using the vincenty distance.

    :param line: Linestring to measure. Coordinates have to be (in the WGS-84 ellipsoid model)
    :return: distance in meters of the linestring.
    """
    coords = line.coords
    if len(coords) < 2:
        return 0.0
    u = coords[0]
    total_distance = 0.0
    for i in range(1, len(coords)):
        v = coords[i]
        total_distance += compute_vincenty(u, v)
        u = v
    return total_distance


def get_new_node_unique_name(graph: nx.Graph, name: str):
    """
    Return a new unique node name from an initial node name. A counter suffix is added at the end if the node name is
    already used.

    :param graph: A given graph
    :param name: A initial node name
    :return: A unique name not in ``graph.nodes()``.
    """
    if name not in graph.nodes():
        return name
    else:
        ct = 1
        unique_name = name
        while unique_name in graph.nodes():
            ct += 1
            unique_name = "%s_%d" % (str(name), ct)
        return unique_name


def euclidian_distance_coordinates(c1: Iterable, c2: Iterable) -> float:
    """Return the euclidian distance between the two sets of coordinates."""
    return math.sqrt(sum(((i - j) ** 2 for i, j in zip(c1, c2))))


def euclidian_distance(p1: Point, p2: Point) -> float:
    """
    Return the euclidian distance between the two points

    :param p1: The first shapely Point
    :param p2: The second shapely Point
    :return: The euclidian distance
    """
    return euclidian_distance_coordinates((p1.x, p1.y), (p2.x, p2.y))


def fill_edges_missing_geometry_attributes(graph: GeoGraph):
    """
    Add a geometry attribute to the edges that don't have any. The created geometry is a straight line between the
    two nodes.

    :param graph: graph to fill
    """
    edges = dict(graph.edges)
    nodes = dict(graph.nodes)
    for s in edges:
        if graph.edges_geometry_key not in edges[s]:
            if graph.nodes_geometry_key in nodes[s[0]] and graph.nodes_geometry_key in nodes[s[1]]:
                p1 = nodes[s[0]][graph.nodes_geometry_key]
                p2 = nodes[s[1]][graph.nodes_geometry_key]
                graph.edges[s][graph.edges_geometry_key] = LineString([p1, p2])


def fill_length_attribute(graph: GeoGraph, attribute_name="length", only_missing=True):
    """
    Fill the ``'length'`` attribute of the given networkX Graph. The length is computed in meters using the vincenty
    formula. Method won't be consistent if the graph crs is not WGS84.

    :param graph: graph to fill
    :param attribute_name: The length attribute name to set
    :param only_missing: Compute the length only if the attribute is missing
    :return: None
    """
    if not crs_equals(graph.crs, settings.USED_CRS):
        raise ValueError("Impossible to compute distance for graph with different"
                         " crs than '%s'. Graph crs: '%s' " % (str(settings.USED_CRS), str(graph.crs)))
    edges_geometry = nx.get_edge_attributes(graph, graph.edges_geometry_key)
    for e in edges_geometry:
        edge_data = graph.edges[e]
        if (not only_missing) or attribute_name not in edge_data:
            edge_data[attribute_name] = measure_line_distance(edges_geometry[e])


def fill_elevation_attribute(graph: GeoGraph, attribute_name="elevation[m]", only_missing=True):
    """
    Fill the ``elevation[m]`` attribute on nodes of the given geograph. The elevation is found with the `srtm` package.
    Graph crs has to be WGS84 standard, otherwise elevation data won't be consistent.

    :param graph: GeoGraph to modify
    :param attribute_name: Attribute to fill
    :param only_missing: Get the elevation and set it only if the node attribute is missing.
    """
    if srtm is None:
        raise ImportError("Impossible to get elevation data, `srtm` package not found.")
    elevation_data = srtm.get_data()
    for n, data in graph.nodes(data=True):
        if (not only_missing) or attribute_name not in data:
            longitude, latitude = tuple(*data[graph.nodes_geometry_key].coords)
            elevation = elevation_data.get_elevation(latitude, longitude)
            if elevation is not None:
                data[attribute_name] = elevation


def join_lines_extremity_to_nodes_coordinates(graph: GeoGraph):
    """
    Modify the edges geometry attribute so that lines extremities match with nodes coordinates.

    :param graph: A geograph to modify
    """
    edges_geometry = nx.get_edge_attributes(graph, graph.edges_geometry_key)
    for e in edges_geometry:
        line = edges_geometry[e]
        to_replace = False
        first_node_coords = graph.get_node_coordinates(e[0])
        if not coordinates_almost_equal(line.coords[0], first_node_coords):
            line = insert_point_in_line(line, first_node_coords, 0)
            to_replace = True
        second_node_coords = graph.get_node_coordinates(e[1])
        if not coordinates_almost_equal(line.coords[-1], second_node_coords):
            line = insert_point_in_line(line, second_node_coords, len(line.coords))
            to_replace = True
        if to_replace:
            graph.edges[e][graph.edges_geometry_key] = line


def get_line_start(graph: GeoGraph, e, line):
    """For a given edge, return the node constituting the line start with a closest node rule."""
    uxy = graph.get_node_coordinates(e[0])
    vxy = graph.get_node_coordinates(e[1])
    first_extremity = line.coords[0]
    last_extremity = line.coords[-1]
    u_e1 = euclidian_distance_coordinates(uxy, first_extremity)
    v_e1 = euclidian_distance_coordinates(vxy, first_extremity)
    u_e2 = euclidian_distance_coordinates(uxy, last_extremity)
    v_e2 = euclidian_distance_coordinates(vxy, last_extremity)
    if u_e1 > v_e1 and u_e2 < v_e2:  # u is closer to e2 and v is closer to e1
        return e[1]
    elif u_e1 < v_e1 and u_e2 < v_e2:  # u is closer to e1 and e2 than v
        if u_e1 > u_e2:
            return e[1]
    elif u_e1 > v_e1 and u_e2 > v_e2:  # v in closer to e1 and e2 than u
        if v_e1 < v_e2:
            return e[1]
    else:
        return e[0]


def get_line_ordered_edge(graph: GeoGraph, e, line):
    """Return the given edge with the first node of the edge representing the first line point and the second node
    the last edge point. The closest node rule is applied."""
    if get_line_start(graph, e, line) != e[0]:
        return (e[1], e[0], *e[2:])
    else:
        return e


def order_well_lines(graph: GeoGraph):
    """
    Try to order well each geometry attribute of edges so that the first coordinates of the line string are the
    coordinates of the first vertex of the edge. The closest node rule is applied. If the graph is not oriented, the
    modification will be inconsistent (nodes declaration in edges views are not ordered). Euclidian distance is used
    here.

    :param graph: Graph on which to apply the ordering step. Modification is inplace.
    :return: None
    """
    line_strings = nx.get_edge_attributes(graph, graph.edges_geometry_key)
    for e, line in line_strings.items():
        if get_line_start(graph, e, line) != e[0]:
            graph.edges[e][graph.edges_geometry_key] = LineString(reversed(line.coords))


def stringify_nodes(graph: nx.Graph, copy=True):
    """Modify the graph node names into strings."""
    nx.relabel_nodes(graph, {n: str(n) for n in graph.nodes}, copy)


def is_nan(val) -> bool:
    return val is np.nan or val != val


def rename_nodes_attribute(graph: nx.Graph, old_name, new_name):
    """Rename nodes attribute defined by its old name to a new name."""
    for n, d in graph.nodes(data=True):
        if old_name in d:
            d[new_name] = d.pop(old_name)


def rename_edges_attribute(graph: nx.Graph, old_name, new_name):
    """Rename edges attribute defined by its old name to a new name."""
    for u, v, d in graph.edges(data=True):
        if old_name in d:
            d[new_name] = d.pop(old_name)


def hard_write_spatial_keys(graph: GeoGraph):
    """Write spatial keys in the graph attribute, so that if the default keys are used, they are propagated for
    special operations (e.g. composing graphs)."""
    for k, v in graph.get_spatial_keys().items():
        setattr(graph, k, v)


def compose(G: GeoGraph, H: GeoGraph) -> GeoGraph:
    """Return a new graph of G composed with H. Makes sure the returned graph is consistent with respect to the spatial
    keys. (See ``nx.compose``). According to the priority rule in networkX, attributes from ``H`` take precedent over
    attributes from G."""
    hard_write_spatial_keys(G)
    hard_write_spatial_keys(H)
    R = nx.compose(G, H)
    if G.nodes_geometry_key != H.nodes_geometry_key:
        rename_nodes_attribute(R, G.nodes_geometry_key, H.nodes_geometry_key)
    if G.edges_geometry_key != H.edges_geometry_key:
        rename_edges_attribute(R, G.edges_geometry_key, H.edges_geometry_key)
    return R


def geographical_distance(graph: GeoGraph, node1, node2, method="vincenty") -> float:
    """Return the geographical distance between the two given nodes.

    :param graph: Geograph
    :param node1: First node label
    :param node2: Second node label
    :param method: "vincenty", "euclidian", "great_circle"
    :return: Distance between nodes, unit depends on the method.
    """
    point1 = graph.nodes[node1][graph.nodes_geometry_key]
    point2 = graph.nodes[node2][graph.nodes_geometry_key]
    if method == "vincenty":
        return compute_vincenty_from_points(point1, point2)
    elif method == "euclidian":
        return euclidian_distance(point1, point2)
    else:
        raise ValueError("Unknown distance method: '%s'" % str(method))


def get_graph_bounding_box(graph: GeoGraph):
    """Return the bounding box coordinates of the given GeoGraph. It takes into account nodes and edges geometries."""
    nodes = list(graph.get_nodes_as_point_series())
    x_min, y_min, x_max, y_max = MultiPoint(nodes).bounds
    bb = [[x_min, y_min], [x_max, y_max]]
    edges = graph.get_edges_as_line_series()
    if len(edges) > 0:
        edges_bounds = edges.bounds
        x_e_min, y_e_min = edges_bounds["minx"].min(), edges_bounds["miny"].min()
        x_e_max, y_e_max = edges_bounds["maxx"].max(), edges_bounds["maxy"].max()
        if x_min > x_e_min:
            bb[0][0] = x_e_min
        if y_min > y_e_min:
            bb[0][1] = y_e_min
        if x_max < x_e_max:
            bb[1][0] = x_e_max
        if y_max < y_e_max:
            bb[1][1] = y_e_max
    return bb
