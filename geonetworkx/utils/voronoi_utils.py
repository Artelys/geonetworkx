"""
    File name: VoronoiParser
    Author: Artelys - Hugo Chareyre
    Date last modified: 21/11/2018
    Python Version: 3.6
"""
import numpy as np
from scipy.spatial import Voronoi, ConvexHull
from shapely.geometry import MultiLineString, LineString, MultiPoint, Point, box, Polygon, MultiPolygon
from shapely.ops import linemerge
import geopandas as gpd
from typing import Union
from collections import defaultdict
import pyvoronoi


GenericLine = Union[LineString, MultiLineString]


class VoronoiParser:
    """Add-on for the scipy.spatial Voronoi tool. It computes the voronoi cells within a bounding box."""

    def __init__(self, points, bounding_box_coords):
        """Constructor of the voronoi parser. It takes the points to compute and a bounding box representing the study
        area."""
        convex_hull = ConvexHull(points)
        self.hull_polygon = Polygon(
            [(x, y) for x, y in zip(points[convex_hull.vertices, 0], points[convex_hull.vertices, 1])])
        self.voronoi_obj = Voronoi(points)
        self.ridges_coords = None
        self.all_regions_coords = None
        self.compute_ridges_coords()
        diagonal_length = np.linalg.norm(np.array(bounding_box_coords[0]) - np.array(bounding_box_coords[1]))
        self.parse_regions(eta=diagonal_length)
        self.bounding_box = box(bounding_box_coords[0][0], bounding_box_coords[0][1],
                                bounding_box_coords[1][0], bounding_box_coords[1][1])

    def compute_ridges_coords(self):
        """Computes the coordinates of the voronoi ridges. It uses the midpoint and the ridge vertex as vector basis. If
        the ridge vertex is not in the convex hull, then we reverse the vector direction so that the cell has the right
        shape."""
        self.ridges_coords = {}
        for p1, p2 in self.voronoi_obj.ridge_dict:
            ridge_vertices = self.voronoi_obj.ridge_dict[(p1, p2)]
            if ridge_vertices[0] < 0 <= ridge_vertices[1]:
                c1 = self.voronoi_obj.vertices[ridge_vertices[1]]
                c2 = (self.voronoi_obj.points[p1] + self.voronoi_obj.points[p2]) / 2.0
            elif ridge_vertices[1] < 0 <= ridge_vertices[0]:
                c1 = self.voronoi_obj.vertices[ridge_vertices[0]]
                c2 = (self.voronoi_obj.points[p1] + self.voronoi_obj.points[p2]) / 2.0
            else:
                continue
            ridge_coords = c2 - c1
            if not self.hull_polygon.intersects(Point(c1)):
                ridge_coords *= -1
            ridge_coords /= np.linalg.norm(ridge_coords)
            self.ridges_coords[(p1, p2)] = ridge_coords

    def parse_regions(self, eta=1.0):
        """Parsing of the voronoi regions using the coordinates of the voronoi ridges. It sets an extremity point at
        infinity (represented by eta) along the ridge. In case of duplicate points, one of the cell will be empty and
        the other won't be."""
        nb_points = len(self.voronoi_obj.points)
        self.all_regions_coords = []
        for p in range(nb_points):
            region = self.voronoi_obj.regions[self.voronoi_obj.point_region[p]]
            region_coords = []
            for ix in range(len(region)):
                v1 = region[ix]
                if v1 >= 0:
                    region_coords.append([self.voronoi_obj.vertices[v1][0], self.voronoi_obj.vertices[v1][1]])
                    continue
                if (ix + 1) < len(region):
                    v2 = region[ix + 1]
                else:
                    v2 = region[0]
                if (ix - 1) >= 0:
                    v0 = region[ix - 1]
                else:
                    v0 = region[-1]
                if v0 != v2:
                    first_ridge_candidates = [(p1, p2) for (p1, p2), vs in self.voronoi_obj.ridge_dict.items()
                                              if (p in (p1, p2)) and (vs == [v0, -1] or vs == [-1, v0])]
                    second_ridge_candidates = [(p1, p2) for (p1, p2), vs in self.voronoi_obj.ridge_dict.items()
                                              if (p in (p1, p2)) and (vs == [v2, -1] or vs == [-1, v2])]
                    if len(first_ridge_candidates) == 0 or len(second_ridge_candidates) == 0:
                        # It means there is a duplicate in points
                        region_coords = []  # This point region will be empty, but the duplicate region won't.
                        break
                    first_ridge = first_ridge_candidates[0]
                    second_ridge = second_ridge_candidates[0]
                else:
                    matching_ridges = [(p1, p2) for p1, p2 in self.voronoi_obj.ridge_dict if (p1 == p or p2 == p) and
                                       (self.voronoi_obj.ridge_dict[(p1, p2)] == [v0, -1] or
                                        self.voronoi_obj.ridge_dict[(p1, p2)] == [-1, v0])]
                    first_ridge = matching_ridges[0]
                    second_ridge = matching_ridges[1]
                first_interpolated_coord = self.voronoi_obj.vertices[v0] + self.ridges_coords[first_ridge] * eta
                second_interpolated_coord = self.voronoi_obj.vertices[v2] + self.ridges_coords[second_ridge] * eta
                region_coords.append([first_interpolated_coord[0], first_interpolated_coord[1]])
                region_coords.append([second_interpolated_coord[0], second_interpolated_coord[1]])
            self.all_regions_coords.append(np.array(region_coords))

    def get_regions_as_polygons(self) -> list:
        """Collection of all the voronoi cells coordinates and creation of shapely polygon with a bounding box trimming
        step."""
        all_polygons = []
        for region in self.all_regions_coords:
            if len(region) > 0:
                polygon = Polygon(region)
                trimmed_polygon = polygon.intersection(self.bounding_box)
                all_polygons.append(trimmed_polygon)
            else:
                all_polygons.append(Polygon())
        return all_polygons

    def get_regions_as_gdf(self, crs=None) -> gpd.GeoDataFrame:
        """Collect all the voronoi cells as shapely polygons and return them as a GeoDataFrame."""
        voronoi_cells_gdf = gpd.GeoDataFrame(columns=['PointId', 'geometry'], crs=crs)
        all_polygons = self.get_regions_as_polygons()
        for point_id, polygon in enumerate(all_polygons):
            voronoi_cells_gdf.loc[len(voronoi_cells_gdf)] = [point_id, polygon]
        return voronoi_cells_gdf




class PyVoronoiHelper:
    """Add-on for the pyvoronoi (boost voronoi) tool. It computes the voronoi cells within a bounding box."""

    def __init__(self, points: list, segments: list, bounding_box_coords: list, scaling_factor=100000.0):
        self.pv = pyvoronoi.Pyvoronoi(scaling_factor)
        for p in points:
            self.pv.AddPoint(p)
        for s in segments:
            self.pv.AddSegment(s)
        points_and_lines = points + [p for l in segments for p in l]
        self.convex_hull = MultiPoint(points_and_lines).convex_hull
        self.convex_hull_centroid = np.array([i for i in self.convex_hull.centroid.coords])
        self.pv.Construct()
        self.discretization_tolerance = 10 / scaling_factor
        self.bounding_box_coords = bounding_box_coords

    def get_cells_as_gdf(self) -> gpd.GeoDataFrame:
        """Returns the voronoi cells in `geodataframe` with a column named `id` referencing the index of the associated
         input geometry."""
        gdf = gpd.GeoDataFrame(columns=["id", "geometry"])
        cells_geometries = self.get_cells_as_polygons()
        gdf["geometry"] = list(cells_geometries.values())
        gdf["id"] = list(cells_geometries.keys())
        return gdf

    def get_cells_as_polygons(self) -> dict:
        """Return the voronoi cells as polygons trimmed with the bounding box."""
        diagonal_length = np.linalg.norm(np.array(self.bounding_box_coords[0]) - np.array(self.bounding_box_coords[1]))
        cells_coordinates = self.get_cells_coordiates(eta=diagonal_length,
                                                      discretization_tolerance=self.discretization_tolerance)
        bounding_box = box(self.bounding_box_coords[0][0], self.bounding_box_coords[0][1],
                           self.bounding_box_coords[1][0], self.bounding_box_coords[1][1])
        cells_as_polygons = dict()
        for i, coords in cells_coordinates.items():
            if len(coords) > 2:
                polygon = Polygon(coords)
                if not polygon.is_valid:
                    polygon = polygon.buffer(0.0)  # TODO: this doesn't fix bowtie case
                trimmed_polygon = polygon.intersection(bounding_box)
                cells_as_polygons[i] = trimmed_polygon
        return cells_as_polygons

    def get_cells_coordiates(self, eta=1.0, discretization_tolerance=0.05) -> dict:
        """"Parse the results of ``pyvoronoi`` to compute the voronoi cells coordinates. The infinite ridges are
        projected at a ``eta`` distance in the ridge direction.

        :param eta: Distance for infinite ridges projection.
        :param discretization_tolerance: Discretization distance for curved edges.
        :return: A dictionary mapping the cells ids and their coordinates.
        """
        vertices = self.pv.GetVertices()
        cells = self.pv.GetCells()
        edges = self.pv.GetEdges()
        cells_coordinates = dict()
        for c in cells:
            cell_coords = []
            for e in c.edges:
                edge = edges[e]
                start_vertex = vertices[edge.start]
                end_vertex = vertices[edge.end]
                if edge.is_linear:
                    if edge.start != -1 and edge.end != -1:
                        self.add_polygon_coordinates(cell_coords, [start_vertex.X, start_vertex.Y])
                        self.add_polygon_coordinates(cell_coords, [end_vertex.X, end_vertex.Y])
                    else:
                        start_is_infinite = edge.start == -1
                        if start_is_infinite:
                            ridge_vertex = end_vertex
                        else:
                            ridge_vertex = start_vertex
                        ridge_point = np.array([ridge_vertex.X, ridge_vertex.Y])
                        twin_cell = cells[edges[edge.twin].cell]
                        if c.site == twin_cell.site:
                            if c.source_category == 3:
                                 segment = np.array(self.pv.RetriveScaledSegment(cells[edge.cell]))
                                 if start_is_infinite:
                                     second_point = segment[1]
                                     ridge_direction = second_point - ridge_point
                                     if np.linalg.norm(ridge_direction) == 0.0:
                                         ridge_direction = - self.get_orthogonal_direction(segment[1] - segment[0])
                                 else:
                                     first_point = segment[0]
                                     ridge_direction = first_point - ridge_point
                                     if np.linalg.norm(ridge_direction) == 0.0:
                                         ridge_direction = - self.get_orthogonal_direction(segment[1] - segment[0])
                            else:
                                first_point = self.pv.RetrieveScaledPoint(c)
                                ridge_direction = first_point - ridge_point
                                if np.linalg.norm(ridge_direction) == 0.0:
                                    segment = np.array(self.pv.RetriveScaledSegment(c))
                                    ridge_direction = - self.get_orthogonal_direction(segment[1] - segment[0])
                                    centroid_direction = first_point - self.convex_hull_centroid
                                    ridge_direction *= np.sign(np.dot(centroid_direction, ridge_direction))
                        else:
                            first_point = self.pv.RetrieveScaledPoint(c)  # TODO : what if first_point == second_point == ridge_point
                            second_point = self.pv.RetrieveScaledPoint(twin_cell)
                            midpoint = np.array([(first_point[0] + second_point[0]) / 2.0,
                                                 (first_point[1] + second_point[1]) / 2.0])
                            ridge_direction = ridge_point - midpoint
                            if self.convex_hull.intersects(Point(ridge_point)):
                                ridge_direction *= -1
                        ridge_direction_norm = np.linalg.norm(ridge_direction)
                        if ridge_direction_norm != 0.0:
                            ridge_direction /= ridge_direction_norm
                        ridge_limit_point = ridge_point + ridge_direction * eta
                        if start_is_infinite:
                            self.add_polygon_coordinates(cell_coords, [ridge_limit_point[0], ridge_limit_point[1]])
                            self.add_polygon_coordinates(cell_coords, [ridge_vertex.X, ridge_vertex.Y])
                        else:
                            self.add_polygon_coordinates(cell_coords, [ridge_vertex.X, ridge_vertex.Y])
                            self.add_polygon_coordinates(cell_coords, [ridge_limit_point[0], ridge_limit_point[1]])
                else:
                    try:
                        coords_to_add = []
                        for p in self.pv.DiscretizeCurvedEdge(e, discretization_tolerance):
                            coords_to_add.append(p)
                        cell_coords.extend(coords_to_add)
                    except pyvoronoi.UnsolvableParabolaEquation:
                        self.add_polygon_coordinates(cell_coords, [start_vertex.X, start_vertex.Y])
                        self.add_polygon_coordinates(cell_coords, [end_vertex.X, end_vertex.Y])
            cells_coordinates[c.cell_identifier] = cell_coords
        return cells_coordinates

    @staticmethod
    def get_orthogonal_direction(dir: np.array) -> np.array:
        """Return the orthogonal direction of the given direction (2D)."""
        if dir[1] == 0.0:
            return np.array([0.0, 1.0])
        return np.array([1.0, - dir[0] / dir[1]])

    @staticmethod
    def add_polygon_coordinates(coordinates: list, point: list):
        """Add given point to given coordinates list if is not the equal to the last coordinates."""
        if coordinates:
            last_point = coordinates[-1]
            if last_point[0] == point[0] and last_point[1] == point[1]:
                return
        coordinates.append(point)


def split_linestring_as_simple_linestrings(line: GenericLine) -> list:
    """Split a linestring if it is not simple (i.e. it crosses itself)."""
    if not line.is_simple:
        mls = line.intersection(line)
        if line.geom_type == 'LineString' and mls.geom_type == 'MultiLineString':
            mls = linemerge(mls)
        return list(mls)
    else:
        return [line]

def split_as_simple_segments(lines: list, tol=1e-6) -> defaultdict:
    """Split a list of lines to simple segments (linestring composed by two points). All returned segments do not
    crosses except at extremities.

    :param lines: List of lines to split
    :param tol: Tolerance to test if a line is a sub line of another one.
    :return: A dictionary mapping for each input line index, the list of simple segments.
    """
    split_lines = defaultdict(list)
    all_split_lines = split_linestring_as_simple_linestrings(MultiLineString(lines))
    lines_stack = [(i, l) for i, l in enumerate(lines)]
    for sub_line in all_split_lines:
        for j, (i, line) in enumerate(lines_stack):
            if line.buffer(tol, 1).contains(sub_line):
                split_lines[i].append(sub_line)
                del lines_stack[j]
                lines_stack.insert(0, (i, line))
                break
    return split_lines


def compute_voronoi_cells_from_lines(lines: list, scaling_factor=1e6) -> gpd.GeoDataFrame:
    """Compute the voronoi cells of given generic lines. Input linestrings can be not simple.

    :param lines: List of ``LineString``
    :param scaling_factor: Resolution for the voronoi cells computation (Two points will be considered equal if their
        coordinates are equal when rounded at ``1/scaling_factor``).
    :return: A `GeoDataFrame` with cells geometries. A column named `id` referencing the index of the associated
        input geometry.
    """
    simple_segments_mapping = split_as_simple_segments(lines, 1 / scaling_factor)
    all_segments = [list(s.coords) for i in range(len(lines)) for s in simple_segments_mapping[i]]
    bounds = MultiLineString(lines).bounds
    bb = [[bounds[0], bounds[1]], [bounds[2], bounds[3]]]
    pvh = PyVoronoiHelper([], segments=all_segments, bounding_box_coords=bb, scaling_factor=scaling_factor)
    gdf = pvh.get_cells_as_gdf()
    gdf = gdf[list(map(lambda i: isinstance(i, Polygon), gdf["geometry"]))].copy()
    gdf["site"] = [pvh.pv.GetCell(c).site for c in gdf["id"]]
    gdf_lines = gpd.GeoDataFrame(columns=["id", "geometry"])
    from shapely.ops import cascaded_union
    from shapely.geometry import GeometryCollection
    ct = 0
    for i, line in enumerate(lines):
        line_polygons = []
        for s in simple_segments_mapping[i]:
            for p in gdf[gdf["site"] == ct]["geometry"]:
                line_polygons.append(p)
            ct += 1
        merged_polygon = cascaded_union(line_polygons)
        if isinstance(merged_polygon, GeometryCollection):
            if len(merged_polygon) == 0:
                continue
            merged_polygon = MultiPolygon(merged_polygon)
        gdf_lines.loc[len(gdf_lines)] = [i, merged_polygon]
    return gdf_lines

