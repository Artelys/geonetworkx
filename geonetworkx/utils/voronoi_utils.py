"""
    File name: VoronoiParser
    Author: Artelys - Hugo Chareyre
    Date last modified: 21/11/2018
    Python Version: 3.6
"""
import numpy as np
from shapely.geometry import MultiLineString, LineString, box, Polygon, MultiPolygon
from shapely.ops import linemerge
import geopandas as gpd
from typing import Union
from collections import defaultdict
import pyvoronoi


GenericLine = Union[LineString, MultiLineString]


class PyVoronoiHelper:
    """Add-on for the pyvoronoi (boost voronoi) tool. It computes the voronoi cells within a bounding box."""

    def __init__(self, points: list, segments: list, bounding_box_coords: list, scaling_factor=100000.0):
        self.pv = pyvoronoi.Pyvoronoi(scaling_factor)
        for p in points:
            self.pv.AddPoint(p)
        for s in segments:
            self.pv.AddSegment(s)
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
                start_vertex = vertices[edge.start] if edge.start != -1 else None
                end_vertex = vertices[edge.end] if edge.end != -1 else None
                if edge.start == -1 or edge.end == -1:
                    self.clip_infinite_edge(cell_coords, edge, eta)
                else:
                    if edge.is_linear:
                        self.add_polygon_coordinates(cell_coords, [start_vertex.X, start_vertex.Y])
                        self.add_polygon_coordinates(cell_coords, [end_vertex.X, end_vertex.Y])
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


    def clip_infinite_edge(self, cell_coords: list, edge: pyvoronoi.Edge, eta: float):
        """Fill infinite edge coordinate by placing the infinite vertex to a ``eta`` distance of the known vertex."""
        cell = self.pv.GetCell(edge.cell)
        twin_edge = self.pv.GetEdge(edge.twin)
        twin_cell = self.pv.GetCell(twin_edge.cell)
        # Infinite edges could not be created by two segment sites.
        if cell.contains_point and twin_cell.contains_point:
            first_point = self.pv.RetrieveScaledPoint(cell)
            second_point = self.pv.RetrieveScaledPoint(twin_cell)
            origin = np.array([(first_point[0] + second_point[0]) / 2.0,
                                 (first_point[1] + second_point[1]) / 2.0])
            ridge_direction = np.array([first_point[1] - second_point[1],
                                        second_point[0] - first_point[0]])
        else:
            if cell.contains_segment:
                origin = np.array(self.pv.RetrieveScaledPoint(twin_cell))
                segment = np.array(self.pv.RetriveScaledSegment(cell))
            else:
                origin = np.array(self.pv.RetrieveScaledPoint(cell))
                segment = np.array(self.pv.RetriveScaledSegment(twin_cell))
            dx = segment[1][0] - segment[0][0]
            dy = segment[1][1] - segment[0][1]
            if (np.linalg.norm(segment[0] - origin) == 0.0) != cell.contains_point:
                ridge_direction = np.array([dy, -dx])
            else:
                ridge_direction = np.array([-dy, dx])
        ridge_direction /= np.linalg.norm(ridge_direction)
        if edge.start == -1:
            ridge_point_projected = origin - ridge_direction * eta
            self.add_polygon_coordinates(cell_coords, ridge_point_projected)
        else:
            start_vertex = self.pv.GetVertex(edge.start)
            self.add_polygon_coordinates(cell_coords, [start_vertex.X, start_vertex.Y])
        if edge.end == -1:
            ridge_point_projected = origin + ridge_direction * eta
            self.add_polygon_coordinates(cell_coords, ridge_point_projected)
        else:
            end_vertex = self.pv.GetVertex(edge.end)
            self.add_polygon_coordinates(cell_coords, [end_vertex.X, end_vertex.Y])


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

def split_as_simple_segments(lines: list, tol=1e-7) -> defaultdict:
    """Split a list of lines to simple segments (linestring composed by two points). All returned segments do not
    cross themselves except at extremities.

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


def compute_voronoi_cells_from_lines(lines: list, scaling_factor=1e7) -> gpd.GeoDataFrame:
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

