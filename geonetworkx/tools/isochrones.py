# -*- coding: utf-8 -*-
"""
    File name: isochrones
    Author: Artelys
    Creation date: 22/05/2019
    Python Version: 3.6
"""
import numpy as np
from geonetworkx.geograph import GeoGraph
from geonetworkx.utils import get_line_start
import geonetworkx as gnx
from shapely.geometry import Polygon, LineString, MultiPolygon
from shapely.ops import cascaded_union
import math
from typing import Union
from scipy.spatial import Delaunay


GenericPolygon = Union[Polygon, MultiPolygon]


def isochrone_polygon(graph: GeoGraph, source, limit, weight="length", tolerance=1e-7) -> GenericPolygon:
    working_graph = graph.copy()
    # Compute the ego-graph
    gnx.add_ego_boundary_nodes(working_graph, source, limit, distance=weight)
    ego_gmg = gnx.extended_ego_graph(working_graph, source, limit, distance=weight)
    # Compute edges voronoi cells
    edge_as_lines = working_graph.get_edges_as_line_series()
    lines = list(edge_as_lines)
    edge_voronoi_cells = gnx.compute_voronoi_cells_from_lines(lines, scaling_factor=1 / tolerance)
    edge_voronoi_cells.set_index("id", inplace=True)
    # Set ego-graph edges cells
    isochrone_polygons = []
    for e, edge in enumerate(edge_as_lines.index):
        if ego_gmg.has_edge(*edge):
            p = edge_voronoi_cells.at[e]
            if any("boundary" in str(n) for n in edge):  # TODO: find a clean way to do this
                boundary_line = edge_as_lines[edge]
                if get_line_start(working_graph, edge, boundary_line) != edge[0]:
                    boundary_line = LineString(reversed(boundary_line.coords))
                edge_buffer_pol = boundary_edge_buffer(boundary_line)
                isochrone_polygons.append(p.intersection(edge_buffer_pol))
            else:
                isochrone_polygons.append(p)
    # Merge as ischrone polygon
    isochrone_polygon = cascaded_union(isochrone_polygons)
    isochrone_polygon = isochrone_polygon.buffer(tolerance)
    return isochrone_polygon






def get_segment_boundary_buffer_polygon(segment_coords: list, radius: float, residual_radius: float) -> Polygon:
    segment_direction = [segment_coords[1][0] - segment_coords[0][0], segment_coords[1][1] - segment_coords[0][1]]
    orthogonal_dir = np.array([- segment_direction[1], segment_direction[0]])
    orthogonal_dir /= np.linalg.norm(orthogonal_dir)
    top_points = [segment_coords[0] + orthogonal_dir * radius, segment_coords[1] + orthogonal_dir * residual_radius]
    bottom_points = [segment_coords[0] - orthogonal_dir * radius, segment_coords[1] - orthogonal_dir * residual_radius]
    return Polygon(top_points + list(reversed(bottom_points)))




def get_point_boundary_buffer_polygon(point_coords: list, radius: float, segment_direction: list, resolution=16) -> Polygon:
    """Returns a half-disk centered on the given point, with the given radius and having the boundary edge orthogonal to
    the given segment direction. See ``boundary_edge_buffer``."""
    # Segment angle with system coordinates
    phi = math.acos(segment_direction[0] / np.linalg.norm(segment_direction))
    if segment_direction[1] < 0:
        phi *= -1.0
    # Discretization angle
    theta = math.pi / float(resolution)
    # Start angle
    angle = phi - math.pi / 2.0
    coords = []
    for i in range(resolution + 1):
        coords.append(point_coords + radius * np.array([math.cos(angle), math.sin(angle)]))
        angle -= theta
    return Polygon(coords)





def boundary_edge_buffer(line: LineString) -> GenericPolygon:
    """ Return the edge buffer polygon on the oriented line. This represented the area where all points are reachable
    starting from the line first extremity and using the closest edge projection rule."""
    radius = line.length
    residual_radius = radius
    boundary_polygons = []
    for i in range(len(line.coords) - 1):
        segment_coords = np.array([line.coords[i], line.coords[i + 1]])
        residual_radius -= gnx.euclidian_distance_coordinates(segment_coords[0], segment_coords[1])
        boundary_polygon = get_segment_boundary_buffer_polygon(segment_coords, radius, residual_radius)
        boundary_polygons.append(boundary_polygon)
        segment_direction = segment_coords[1] - segment_coords[0]
        boundary_polygons.append(get_point_boundary_buffer_polygon(segment_coords[0], radius, segment_direction))
        radius = residual_radius
    return cascaded_union(boundary_polygons)



def get_alpha_shape_polygon(points: list, quantile: int) -> GenericPolygon:
    """Return the alpha-shape polygon formed by the given points. Alpha parameter is determined using a quantile of
    circumradius of Delaunay triangles.

    :param points: List of input points (2D)
    :param quantile: Quantile on circumradius to determine alpha (100 returns the convex hull,
        0 returns an empty polygon). ``0 <= quantile <= 100``.
    :return: The polygon formed by all triangles having a circumradius inferior or equal to :math:`1/\\alpha`.

    Note that this does not return the exhaustive alpha-shape for low quantiles, the minimum spanning tree LineString
    should be added to the returned polygon.
    """
    tri = Delaunay(np.array(points))
    polygons = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    circum_radius = []
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Lengths of sides of triangle
        a = gnx.euclidian_distance_coordinates(pa, pb)
        b = gnx.euclidian_distance_coordinates(pb, pc)
        c = gnx.euclidian_distance_coordinates(pc, pa)
        # Semiperimeter of triangle
        s = (a + b + c) / 2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        circum_radius.append(circum_r)
        polygons.append(Polygon([pa, pb, pc]))
    inv_alpha = np.percentile(circum_radius, quantile)
    filtered_polygons = [p for i, p in enumerate(polygons) if circum_radius[i] <= inv_alpha]
    return cascaded_union(filtered_polygons)
