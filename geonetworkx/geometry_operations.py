# -*- coding: utf-8 -*-
import numpy as np
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString
from scipy.spatial import cKDTree
from collections import defaultdict
from typing import Union, Iterable
import pyproj
import geonetworkx.settings as settings
import geonetworkx as gnx


PointCoordinatesLike = Iterable[float]
PointsCoordinatesLike = Union[Iterable[PointCoordinatesLike], MultiPoint]


class Extremity:
    """Represents an extremity of a line. It's useful to parse and deal with lines given as input."""
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
    """Merge two lines (line1 and line2) with the given extremities (e1 and e2)."""
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


def merge_two_lines_with_closest_extremities(first_line: LineString, second_line: LineString) -> LineString:
    """Merge two lines with their closest extremities. Euclidian distance is used here."""
    combinations = [(0, 0), (-1, 0), (0, -1), (-1, -1)]
    extremities = np.array([[first_line.coords[i1], second_line.coords[i2]] for (i1, i2) in combinations])
    distances = np.linalg.norm(extremities[:, 0, :] - extremities[:, 1, :], axis=1)
    merge_couple_index = int(np.argmin(distances))
    merge_couple = combinations[merge_couple_index]
    if merge_couple == (-1, 0):
        merged_line = LineString(list(first_line.coords) + list(second_line.coords))
    elif merge_couple == (0, 0):
        merged_line = LineString(list(reversed(first_line.coords)) + list(second_line.coords))
    elif merge_couple == (0, -1):
        merged_line = LineString(list(reversed(first_line.coords)) + list(reversed(second_line.coords)))
    else:  # merge_couple == (-1, -1)
        merged_line = LineString(list(first_line.coords) + list(reversed(second_line.coords)))
    return merged_line


def get_shape_extremities(shape: LineString, shape_id: int):
    """Return the extremities of a shape in the network_shapes_gdf."""
    assert isinstance(shape, LineString), "Shape must be a LineString"
    first_vertex = shape.coords[0]
    last_vertex = shape.coords[-1]
    e1 = Extremity(shape_id, 0, first_vertex)
    e2 = Extremity(shape_id, -1, last_vertex)
    e1.related_extremity = e2
    e2.related_extremity = e1
    return e1, e2


def convert_multilinestring_to_linestring(gdf: gpd.GeoDataFrame) -> int:
    """Convert all geometry attribute being a 'MultiLineString' to a 'LineString'. The created line is a merge of all
     sublines.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        A GeoDataFrame with a 'geometry' column to modify

    Returns
    -------
    int
        The number of converted 'MultiLineString'

    Raises
    ------
    RuntimeError
        If an input shape is not a LineString or a MultiLineString

    """
    nb_converted_multilinestring = 0
    for s in gdf.index:
        if isinstance(gdf.at[s, 'geometry'], LineString):
            continue
        elif isinstance(gdf.at[s, 'geometry'], MultiLineString):
            all_lines = gdf.at[s, 'geometry'].geoms
            new_line = all_lines[0]
            for i in range(1, len(all_lines)):
                new_line = merge_two_lines_with_closest_extremities(new_line, all_lines[i])
            gdf.at[s, 'geometry'] = new_line
            nb_converted_multilinestring += 1
        else:
            raise RuntimeError("Unknown shape type for shape input at index : ", s)
    return nb_converted_multilinestring


def discretize_line(line: LineString, discretization_tol) -> list:
    """Takes a shapely LineString and discretize it into a list of shapely Points. Each point is at most at the
    discretization tolerance distance of the following point.

    Parameters
    ----------
    line : LineString
        Line to discretize
    discretization_tol : float
        Maximum distance between two points on the line.

    Returns
    -------
    list
        An ordered list of shapely Point

    See Also
    --------
    discretize_lines

    """
    if discretization_tol <= 0.0:
        raise ValueError("Discretization tolerance must be strictly positive.")
    points_list = []
    current_dist = discretization_tol
    line_length = line.length
    points_list.append(Point(list(line.coords)[0]))
    while current_dist < line_length:
        points_list.append(line.interpolate(current_dist))
        current_dist += discretization_tol
    points_list.append(Point(list(line.coords)[-1]))
    return points_list


def discretize_lines(lines: Iterable[LineString], discretization_tol):
    """Discretize some line into points.

    Parameters
    ----------
    lines: Iterable[LineString] :
        Lines to discretize
    discretization_tol : float
        Maximum distance between two points on the line.

    Returns
    -------
    MultiPoint and defaultdict
        Return all the discretized points as a shapely MultiPoint and a dictionary to map the discretized points for
        each line.

    See Also
    --------
    discretize_line

    """
    points_line_association = defaultdict(list)
    all_points = []
    nb_points = 0
    for line_index, line in enumerate(lines):
        discretized_points = discretize_line(line, discretization_tol)
        nb_discretized_points = len(discretized_points)
        all_points.extend(discretized_points)
        points_line_association[line_index] = list(range(nb_points, nb_points + nb_discretized_points))
        nb_points += nb_discretized_points
    all_points_as_mp = MultiPoint(all_points)
    return all_points_as_mp, points_line_association


def get_closest_point_from_points(points_from: PointsCoordinatesLike, points_to: list = None, kd_tree: cKDTree = None):
    """Compute the closest point among the ``points_from`` list for each point in the ``points_to`` list.

    Parameters
    ----------
    points_from : PointsCoordinatesLike
        Iterable of points coordinates
    points_to : list or None
        Iterable of points coordinates (Default value = None)
    kd_tree : cKDTree
        a constructed kd tree representing ``points_from`` (Default value = None)

    Returns
    -------
    array of floats
        distances
    ndarray of ints
        indexes

    """
    if points_to is None and kd_tree is None:
        raise ValueError("Must provide at least argument 'points_to' or 'kd_tree'")
    if kd_tree is None:
        kd_tree = cKDTree(points_to)
    return kd_tree.query(points_from)


def get_closest_point_from_line(line_from: LineString, discretization_tol: float,
                                points_to: list = None, kd_tree: cKDTree = None):
    """Return the closest point from a given line and its distance.

    Parameters
    ----------
    line_from : LineString
        A shapely LineString (Default value = None)
    discretization_tol : float
        Maximum distance between two discretized points on the line.
    points_to : list
        A list of points among which the closest to the line has to be found (optional is ``kdtree`` is
        given)
    kd_tree : cKDTree
        A kd-tree representing the points among which the closest to the line has to be found (optional if
        ``points_to`` is given) (Default value = None)

    Returns
    -------
    float
        closest distance
    int
        index of the closest point

    """
    if points_to is None and kd_tree is None:
        raise ValueError("Must provide at least argument 'points_to' or 'kd_tree'")
    if kd_tree is None:
        kd_tree = cKDTree(points_to)
    discretized_points = discretize_line(line_from, discretization_tol)
    distances, closest_points_indexes = kd_tree.query(MultiPoint(discretized_points))
    smallest_distance_index = np.argmin(distances)
    return distances[smallest_distance_index], closest_points_indexes[smallest_distance_index]


def get_closest_point_from_multi_shape(multi_shape, points_to=None, kd_tree=None):
    """Computes the closest point to the multi shape (i.e. the point that has the smallest projection distance on the
    entire multi shape object.

    Parameters
    ----------
    multi_shape :  MultiPoint or MultiLineString
        The multi shape object can be any shapely object among: MultiPoint, MultiLineString
    points_to : list
        A list of points among which to find the closest to the multi shape (Default value = None)
    kd_tree : cKDTree
        A kdtree representing the points among which the closest to the multishape has to be found (optional
        if 'points_to' is given) (Default value = None)

    Returns
    -------
    float
        distance
    int
        index of the closest point

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
                                 kd_tree: cKDTree = None):
    """Compute the closest point to the given shape.

    Parameters
    ----------
    shape : Point, MultiPoint, LineString or MultiLineString
        Any shapely shape
    points_to : list
        A list of points among which to find the closest to the multi shape (Default value = None)
    kd_tree : cKDTree
        A kdtree representing the points among which the closest to the shape has to be found (optional if
        'points_to' is given) (Default value = None)

    Returns
    -------
    float
        distance
    int
        index of the closest point

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
    """Compute the closest point for each given shape.

    Parameters
    ----------
    shapes_from :
        An iterable of shapes (Point, MultiPoint, LineString, MultiLineString)
    points_to : list
        A list of points among which to find the closest to the multi shape

    Returns
    -------
    float
        distance
    int
        index of the closest point

    """
    results = []
    kd_tree = cKDTree(points_to)
    for shape in shapes_from:
        distance, closest_points_index = get_closest_point_from_shape(shape, kd_tree=kd_tree)
        results.append((distance, closest_points_index))
    return results


def get_closest_line_from_point(point_from: PointCoordinatesLike,
                                lines_to=None,
                                discretization_tol=None,
                                kd_tree=None,
                                points_line_association=None):
    """Find the closest line from a given point.

    Parameters
    ----------
    point_from : PointCoordinatesLike
        Point coordinate to find the closest line.
    lines_to : list
        Group of lines among which the closest has to be found (optional if ``kdtree`` and
        ``points_line_association`` are given). (Default value = None)
    discretization_tol: float
        Maximum distance between discretized points (optional if ``kdtree`` and
        ``points_line_association`` are given). (Default value = None)
    kd_tree : cKDTree
        An optional pre-computed kd_tree of discretized lines. (Default value = None)
    points_line_association : dict
        An optional pre-computed dictionary matching lines and discretized points. (Default value = None)

    Returns
    -------
    float
        distance
    int
        index of the closest line

    """
    if lines_to is None and kd_tree is None:
        raise ValueError("Must provide at least argument 'points_to' or 'kd_tree'")
    elif kd_tree is None and points_line_association is None:
        raise ValueError("If a kd-tree is given, a point line association dictionary must provided")
    if kd_tree is None:
        if discretization_tol is None:
            raise ValueError("If no kd-tree is given, a discretization tolerance must be provided.")
        points_to, points_line_association = discretize_lines(lines_to, discretization_tol)
        kd_tree = cKDTree(points_to)
    distance, closest_point_index = get_closest_point_from_points([point_from], kd_tree=kd_tree)
    for line_index in points_line_association:
        if closest_point_index in points_line_association[line_index]:
            return distance, line_index
    assert False


def get_closest_line_from_points(points_from, lines_to, discretization_tol):
    """Find the closest line for each given points.

    Parameters
    ----------
    points_from : list
        Points coordinates.
    lines_to : list
        Group of lines among which the closest has to be found.
    discretization_tol: float
        Maximum distance between discretized points

    Returns
    -------
    list
        A list of closest lines indexes.

    """
    points_to, points_line_association = discretize_lines(lines_to, discretization_tol)
    kd_tree = cKDTree(points_to)
    lines_indexes = []
    for point in points_from:
        result = get_closest_line_from_point(point, kd_tree=kd_tree, points_line_association=points_line_association)
        lines_indexes.append(result[1])
    return lines_indexes


def get_polygons_neighborhood(polygons):
    """Returns for each polygon a set of intersecting polygons."""
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
    """Cuts a line in two at a distance from its starting point."""
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
    """Returns true if the two given list of coordinates equals within a given tolerance.

    Parameters
    ----------
    c1: Iterable
        First point coordinates
        
    c2: Iterable
        Second point coordinates
        
    tolerance : float
         Tolerance comparison (Default value = 1e-8)

    Returns
    -------
    bool
        True if the coordinates almost equal, false otherwise.

    """
    for i, j in zip(c1, c2):
        if abs(i - j) > tolerance:
            return False
    return True


def almost_equally_located(p1: Point, p2: Point, tolerance=1e-8) -> bool:
    """Test if two point are loacated at the same place within a tolerance.

    Parameters
    ----------
    p1: Point
        First point to compare
    p2: Point
        Second point to compare
    tolerance :
         Comparison tolerance (Default value = 1e-8)

    Returns
    -------
    bool
        True if the two points have the same coordinates.
    """
    return coordinates_almost_equal([p1.x, p1.y], [p2.x, p2.y], tolerance)


def insert_point_in_line(line: LineString, point_coords: list, position: int) -> LineString:
    """Insert a new point in a line given its coordinates."""
    new_line_coordinates = line.coords[0:position]
    new_line_coordinates.append((point_coords[0], point_coords[1]))
    new_line_coordinates.extend(line.coords[position:])
    return LineString(new_line_coordinates)


def get_default_discretization_tolerance(crs):
    """Return a discretization tolerance with the right order of magnitude for
    the given crs.

    Examples
    --------
    >>> import geonetworkx as gnx
    >>> print(gnx.get_default_discretization_tolerance("epsg:3857"))
    3.0
    """
    if not gnx.is_null_crs(crs):
        pp_crs = pyproj.CRS(crs)
        if gnx.crs_equals(pp_crs, gnx.WGS84_CRS):
            return 1e-4  # in degree
        elif pp_crs.axis_info is not None and len(pp_crs.axis_info) > 0:
            first_axis = pp_crs.axis_info[0]
            if first_axis.unit_name == "metre":
                return 3.0  # in meters
    raise ValueError("Impossible to provide a valid discretization tolerance"
                     "for the given crs. A custom one can be set.")