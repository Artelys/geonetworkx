import numpy as np
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString
from scipy.spatial import KDTree
from collections import defaultdict
from typing import Union, Iterable
import geonetworkx.settings as settings


PointCoordinatesLike = Iterable[float]
PointsCoordinatesLike = Union[Iterable[PointCoordinatesLike], MultiPoint]


class Extremity:
    """
    Represents an extremity of a line. It's useful to parse and deal with lines given as input.
    """
    nb_extremity = 0

    def __init__(self, shape_id, position, coords):
        self.id = Extremity.nb_extremity
        Extremity.nb_extremity += 1
        self.shape_id = shape_id
        self.position = position
        self.coords = coords
        self.related_extremity = None
        self.close_extremities = []
        self.closest_station = None
        self.is_merged = False
        self.is_added_as_shape = False
        self.matching_items = dict()


def merge_two_shape(e1: Extremity, e2: Extremity, line1: LineString, line2: LineString) -> LineString:
    """
    Merge two lines (line1 and line2) with the given extremities (e1 and e2).

    :param e1: line1 extremity
    :param e2: line2 extremity
    :param line1: first line (shapely LineString)
    :param line2: second line (shapely LineString)
    :return: the resulting merged line (shapely LineString)
    """
    if e1.position == -1 and e2.position == 0:
        new_line = LineString(list(line1.coords) + list(line2.coords))
    elif e1.position == 0 and e2.position == 0:
        new_line = LineString(list(reversed(line1.coords)) + list(line2.coords))
    elif e1.position == 0 and e2.position == -1:
        new_line = LineString(list(reversed(line1.coords)) + list(reversed(line2.coords)))
    elif e1.position == -1 and e2.position == -1:
        new_line = LineString(list(line1.coords) + list(reversed(line2.coords)))
    else:
        assert False, 'Invalid position attribute'
    return new_line


def merge_two_shapes_with_closest_extremities(shape1: LineString, shape2: LineString) -> LineString:
    """
    Merge two lines with their closest extremities

    :param shape1: first line
    :param shape2: second line
    :return: merged shape
    """
    combinations = [(0, 0), (-1, 0), (0, -1), (-1, -1)]
    extremities = [(np.array(shape1.coords[i1]), np.array(shape2.coords[i2])) for (i1, i2) in combinations]
    distances = [np.linalg.norm(e1 - e2) for (e1, e2) in extremities]
    merge_couple = int(np.argmin(distances))
    e1 = Extremity(None, combinations[merge_couple][0], shape1.coords[combinations[merge_couple][0]])
    e2 = Extremity(None, combinations[merge_couple][1], shape1.coords[combinations[merge_couple][1]])
    return merge_two_shape(e1, e2, shape1, shape2)


def merge_edges_connected_with_two_degree_node(graph: nx.Graph, filter=None) -> int:
    """
    Merge all nodes with one incoming edge and one outgoing edge.

    :param graph: The graph to modify
    :param filter: A lambda function indicating if a given node is has to be potentially merge
    :return: The number of merged nodes
    """
    if filter is None:
        two_degree_nodes = [x for x in graph.nodes() if graph.degree(x) == 2]
    else:
        two_degree_nodes = [x for x in graph.nodes() if graph.degree(x) == 2 and filter(x)]
    nb_two_degree_nodes_merged = 0
    for n in two_degree_nodes:
        edges = [e for e in graph.edges if n in [e[0], e[1]]]
        # Skip if edges are parallel edges for MultiGraph
        if (edges[0][0], edges[0][1]) == (edges[1][0], edges[1][1]):
            continue
        # merge with closest extremities
        shape1 = graph.edges[edges[0]]['geometry']
        shape2 = graph.edges[edges[1]]['geometry']
        new_line = merge_two_shapes_with_closest_extremities(shape1, shape2)
        # Remove node from graph and add one replacing edge
        graph.remove_node(n)
        u = edges[0][0] if edges[0][0] != n else edges[0][1]
        v = edges[1][0] if edges[1][0] != n else edges[1][1]
        graph.add_edge(u, v, geometry=new_line)
        nb_two_degree_nodes_merged += 1
    return nb_two_degree_nodes_merged


def get_shape_extremities(shape: LineString, shape_id: int):
    """
    Return the extremities of a shape in the network_shapes_gdf.

    :param shape: LineString on which to parse Extremity objects
    :param shape_id: id of the shape
    :return: Two extremities of the shape.
    """
    assert isinstance(shape, LineString), "Shape must be a LineString"
    first_vertex = shape.coords[0]
    last_vertex = shape.coords[-1]
    e1 = Extremity(shape_id, 0, first_vertex)
    e2 = Extremity(shape_id, -1, last_vertex)
    e1.related_extremity = e2
    e2.related_extremity = e1
    return e1, e2


def convert_multilinestring_to_linestring(gdf: gpd.GeoDataFrame) -> int:
    """
    Convert all geometry attribute being a 'MultiLineString' to a 'LineString'. The created line is a merge of all sub
    lines.

    :param gdf: A GeoDataFrame with a 'geometry' column to modify
    :return: The number of converted 'MultiLineString'
    """
    nb_converted_multilinestring = 0
    for s in gdf.index:
        if isinstance(gdf.at[s, 'geometry'], LineString):
            continue
        elif isinstance(gdf.at[s, 'geometry'], MultiLineString):
            all_lines = gdf.at[s, 'geometry'].geoms
            new_line = all_lines[0]
            for i in range(1, len(all_lines)):
                new_line = merge_two_shapes_with_closest_extremities(new_line, all_lines[i])
            gdf.at[s, 'geometry'] = new_line
            nb_converted_multilinestring += 1
        else:
            raise RuntimeError("Unknown shape type for shape input at index : ", s)
    return nb_converted_multilinestring


def discretize_line(line: LineString):
    """
    Takes a shapely LineString and discretize it into a list of shapely Points. Each point is at most at the
    discretization tolerance distance of the following point.

    :param line: Line to discretize
    :return: An ordered list of shapely Point
    """
    points_list = []
    current_dist = settings.DISCRETIZATION_TOLERANCE
    line_length = line.length
    points_list.append(Point(list(line.coords)[0]))
    while current_dist < line_length:
        points_list.append(line.interpolate(current_dist))
        current_dist += settings.DISCRETIZATION_TOLERANCE
    points_list.append(Point(list(line.coords)[-1]))
    return points_list


def discretize_lines(lines: Iterable[LineString]):
    points_line_association = defaultdict(list)
    all_points = []
    nb_points = 0
    for line_index, line in enumerate(lines):
        discretized_points = discretize_line(line)
        nb_discretized_points = len(discretized_points)
        all_points.extend(discretized_points)
        points_line_association[line_index] = list(range(nb_points, nb_points + nb_discretized_points))
        nb_points += nb_discretized_points
    all_points_as_mp = MultiPoint(all_points)
    return all_points_as_mp, points_line_association


def get_closest_point_from_points(points_from: PointsCoordinatesLike, points_to: list = None, kd_tree: KDTree = None):
    """
    Compute the closest point among the 'points_from' list for each point in the 'points_from' list.

    :param points_from: Iterable of points coordinates
    :param points_to: Iterable of points coordinates
    :param kd_tree: a constructed kd tree representing `points_from`
    :return: tuple: (distances, indexes)
    """
    if points_to is None and kd_tree is None:
        raise ValueError("Must provide at least argument 'points_to' or 'kd_tree'")
    if kd_tree is None:
        kd_tree = KDTree(points_to)
    return kd_tree.query(points_from)


def get_closest_point_from_line(line_from: LineString, points_to: list = None, kd_tree: KDTree = None):
    """
    Return the closest point from a given line and its distance.

    :param line_from: A shapely LineString
    :param points_to: A list of points among which the closest to the line has to be found (optional is `kdtree` is
        given)
    :param kd_tree: A kd-tree representing the points among which the closest to the line has to be found (optional if
        'points_to' is given)
    :return: a couple containing the closest distance and the index of the closest point
    """
    if points_to is None and kd_tree is None:
        raise ValueError("Must provide at least argument 'points_to' or 'kd_tree'")
    if kd_tree is None:
        kd_tree = KDTree(points_to)
    discretized_points = discretize_line(line_from)
    distances, closest_points_indexes = kd_tree.query(MultiPoint(discretized_points))
    smallest_distance_index = np.argmin(distances)
    return distances[smallest_distance_index], closest_points_indexes[smallest_distance_index]


def get_closest_point_from_multi_shape(multi_shape, points_to=None, kd_tree=None):
    """
    Computes the closest point to the multi shape (i.e. the point that has the smallest projection distance on the
    entire multi shape object.

    :param multi_shape: The multi shape object can be any shapely object among: MultiPoint, MultiLineString
    :param points_to:  A list of points among which to find the closest to the multi shape
    :param kd_tree: A kdtree representing the points among which the closest to the multishape has to be found (optional
        if 'points_to' is given)
    :return: A couple containing the distance and the index of the closest point
    """
    if isinstance(multi_shape, MultiPoint):
        distances_and_points_indexes = get_closest_point_from_points(multi_shape, points_to, kd_tree)
    elif isinstance(multi_shape, MultiLineString):
        distances_and_points_indexes = []
        for line in multi_shape:
            distances_and_points_indexes.append(get_closest_point_from_line(line, points_to, kd_tree))
    else:
        raise TypeError("Method not implemented for shape of type '%s', expected MultiPoint"
                        "or MultiLineString" % str(type(multi_shape)))
    smallest_distance_index = np.argmin(d_ix[0] for d_ix in distances_and_points_indexes)
    return distances_and_points_indexes[smallest_distance_index]


def get_closest_point_from_shape(shape: Union[Point, LineString, MultiPoint, MultiLineString],
                                 points_to: Union[MultiPoint, np.ndarray, list] = None,
                                 kd_tree: KDTree = None):
    """
    Compute the closest point to the given shape.

    :param shape: Any shapely shape (Point, MultiPoint, LineString, MultiLineString)
    :param points_to:  A list of points among which to find the closest to the multi shape
    :param kd_tree: A kdtree representing the points among which the closest to the shape has to be found (optional if
        'points_to' is given)
    :return: A couple containing the distance and the index of the closest point
    """
    if isinstance(shape, Point):
        result = get_closest_point_from_points(MultiPoint([shape]), points_to, kd_tree=kd_tree)
        return result[0][0], result[1][0]
    elif isinstance(shape, LineString):
        return get_closest_point_from_line(shape, points_to, kd_tree=kd_tree)
    elif isinstance(shape, (MultiPoint, MultiLineString)):
        return get_closest_point_from_multi_shape(shape, points_to, kd_tree=kd_tree)
    else:
        raise TypeError("Method not implemented for shape of type '%s', expected Point, "
                        "LineString or MultiLineString" % str(type(shape)))


def get_closest_point_from_shapes(shapes_from, points_to):
    """
    Compute the closest point for each given shape.

    :param shapes_from: An iterable of shapes (Point, MultiPoint, LineString, MultiLineString)
    :param points_to:  A list of points among which to find the closest to the multi shape
    :return: A list of couple containing the distance and the index of the closest point
    """
    results = []
    kd_tree = KDTree(points_to)
    for shape in shapes_from:
        distance, closest_points_index = get_closest_point_from_shape(shape, kd_tree=kd_tree)
        results.append((distance, closest_points_index))
    return results


def get_closest_line_from_point(point_from: PointCoordinatesLike,
                                lines_to=None,
                                kd_tree=None,
                                points_line_association=None):
    """
    Find the closest line from a given point.

    :param point_from: Point coordinate to find the closest line.
    :param lines_to: Group of lines among which the closest has to be found (optional if `kdtree` and
                    `points_line_association` are given).
    :param kd_tree: An optional pre-computed kd_tree of discretized lines.
    :param points_line_association: An optional pre-computed dictionary matching lines and discretized points.
    :return: The closest distance and closest line index.
    """
    if lines_to is None and kd_tree is None:
        raise ValueError("Must provide at least argument 'points_to' or 'kd_tree'")
    elif kd_tree is None and points_line_association is None:
        raise ValueError("If a kd-tree is given, a point line association dictionary must provided")
    if kd_tree is None:
        points_to, points_line_association = discretize_lines(lines_to)
        kd_tree = KDTree(points_to)
    distance, closest_point_index = get_closest_point_from_points([point_from], kd_tree=kd_tree)
    for line_index in points_line_association:
        if closest_point_index in points_line_association[line_index]:
            return distance, line_index
    assert False


def get_closest_line_from_points(points_from, lines_to):
    """
    Find the closest line for each given points.

    :param points_from: Points coordinates.
    :param lines_to: Group of lines among which the closest has to be found.
    :return: A list of closest lines indexes.
    """
    points_to, points_line_association = discretize_lines(lines_to)
    kd_tree = KDTree(points_to)
    lines_indexes = []
    for point in points_from:
        result = get_closest_line_from_point(point, kd_tree=kd_tree, points_line_association=points_line_association)
        lines_indexes.append(result[1])
    return lines_indexes


def get_polygons_neighborhood(polygons):
    nb_polygons = len(polygons)
    neighborhoods = [set() for c in range(nb_polygons)]
    for ix1, p1 in enumerate(polygons):
        for ix2 in range(ix1 + 1, nb_polygons):
            p2 = polygons[ix2]
            if p1.intersects(p2):
                neighborhoods[ix1].add(ix2)
                neighborhoods[ix2].add(ix1)
    return neighborhoods


def split_line(line, distance):
    """Cuts a line in two at a distance from its starting point"""
    coords = list(line.coords)
    if distance <= 0.0:
        return [LineString([coords[0]]*2), LineString(line)]
    if distance >= line.length:
        return [LineString(line), LineString([coords[-1]]*2)]
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [LineString(coords[:i + 1]), LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [LineString(coords[:i] + [(cp.x, cp.y)]), LineString([(cp.x, cp.y)] + coords[i:])]


def coordinates_almost_equal(c1: Iterable, c2: Iterable, tolerance=1e-8) -> bool:
    """Return true if the two given set of coordinates are almost equals within a given tolerance"""
    for i, j in zip(c1, c2):
        if abs(i - j) > tolerance:
            return False
    return True


def almost_equally_located(p1: Point, p2: Point, tolerance=1e-8) -> bool:
    """Return True if the two given point are equally located within a given tolerance"""
    return coordinates_almost_equal([p1.x, p1.y], [p2.x, p2.y], tolerance)


def insert_point_in_line(line: LineString, point_coords: list, position: int):
    """Insert a new point in a line given its coordinates"""
    new_line_coordinates = line.coords[0:position]
    new_line_coordinates.append((point_coords[0], point_coords[1]))
    new_line_coordinates.extend(line.coords[position:])
    return LineString(new_line_coordinates)

