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
import geopandas as gpd


GenericPolygon = Union[Polygon, MultiPolygon]


def get_edges_voronoi_cells(graph: GeoGraph, tolerance=1e-7) -> gpd.GeoSeries:
    """Return edge voronoi cells as `GeoSeries`."""
    edge_as_lines = graph.get_edges_as_line_series()
    lines = list(edge_as_lines)
    edge_cells = gnx.compute_voronoi_cells_from_lines(lines, tolerance)
    edge_cells_as_series = gpd.GeoSeries({e: cell for e, cell in zip(edge_as_lines.index, edge_cells)})
    edge_cells_as_series.crs = graph.crs
    return edge_cells_as_series


def get_segment_boundary_buffer_polygon(segment_coords: Union[list, np.array],
                                        radius: float, residual_radius: float) -> Polygon:
    """Return a segment boundary polygon using given radius. It represents all reachable points from the first
    extremity of the segment. The returned polygon is a trapeze. See ``boundary_edge_buffer``."""
    segment_direction = [segment_coords[1][0] - segment_coords[0][0], segment_coords[1][1] - segment_coords[0][1]]
    orthogonal_dir = np.array([- segment_direction[1], segment_direction[0]])
    orthogonal_dir /= np.linalg.norm(orthogonal_dir)
    top_points = [segment_coords[0] + orthogonal_dir * radius, segment_coords[1] + orthogonal_dir * residual_radius]
    bottom_points = [segment_coords[0] - orthogonal_dir * radius, segment_coords[1] - orthogonal_dir * residual_radius]
    return Polygon(top_points + list(reversed(bottom_points)))


def get_point_boundary_buffer_polygon(point_coords: list, radius: float,
                                      segment_direction: list, resolution=16) -> Polygon:
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
    """Return the edge buffer polygon on the oriented line. This represented the area where all points are reachable
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


def isochrone_polygon(graph: GeoGraph, source, limit, weight="length", tolerance=1e-7) -> GenericPolygon:
    """Return a polygon approximating the isochrone set in the geograph.

    :param graph: Graph representing possible routes.
    :param source: Source node from where distance is computed
    :param limit: Isochrone limit (e.g. 100 meters, 5 minutes, depending on ``weight`` unit).
    :param weight: Weight attribute on edges to compute distances (edge weights should be non-negative).
    :param tolerance: Tolerance to compute Voronoi cells.
    :return: A polygon representing all reachable points within the given limit from the source node.
    """
    working_graph = graph.copy()
    # Compute the ego-graph
    gnx.add_ego_boundary_nodes(working_graph, source, limit, distance=weight)
    ego_graph = gnx.extended_ego_graph(working_graph, source, limit, distance=weight)
    # Compute edges voronoi cells
    edge_voronoi_cells = get_edges_voronoi_cells(working_graph, tolerance)
    # Set ego-graph edges cells
    isochrone_polygons = []
    for e, edge in enumerate(edge_voronoi_cells.index):
        if ego_graph.has_edge(*edge):
            p = edge_voronoi_cells.at[e]
            if not graph.has_edge(*edge):  # For boundary edges
                boundary_line = working_graph.edges[edge][working_graph.edges_geometry_key]
                if get_line_start(working_graph, edge, boundary_line) != edge[0]:
                    boundary_line = LineString(reversed(boundary_line.coords))
                edge_buffer_pol = boundary_edge_buffer(boundary_line)
                isochrone_polygons.append(p.intersection(edge_buffer_pol))
            else:
                isochrone_polygons.append(p)
    # Merge as isochrone polygon
    final_polygon = cascaded_union(isochrone_polygons)
    final_polygon = final_polygon.buffer(tolerance)
    return final_polygon


def get_alpha_shape_polygon(points: list, quantile: float) -> GenericPolygon:
    """Return the alpha-shape polygon formed by the given points. Alpha parameter is determined using a quantile of
    circumradius of Delaunay triangles.

    :param points: List of input points (2D)
    :param quantile: Quantile on circumradius to determine alpha (100 returns the convex hull,
        0 returns an empty polygon). ``0 <= quantile <= 100``.
    :return: The polygon formed by all triangles having a circumradius inferior or equal to :math:`1/\\alpha`.

    Note that this does not return the exhaustive alpha-shape for low quantiles, the minimum spanning tree LineString
    should be added to the returned polygon.
    This is adapted from `Sean Gillies code <https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html>`_.
    """
    points = np.asarray(points)
    tri = Delaunay(points)
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
        area = math.sqrt(np.clip(s * (s - a) * (s - b) * (s - c), 0.0, float("inf")))
        if area <= 0.0:
            continue
        circum_r = a * b * c / (4.0 * area)
        circum_radius.append(circum_r)
        polygons.append(Polygon([pa, pb, pc]))
    inv_alpha = np.percentile(circum_radius, quantile)
    filtered_polygons = [p for i, p in enumerate(polygons) if circum_radius[i] <= inv_alpha]
    return cascaded_union(filtered_polygons)


def isochrone_polygon_with_alpha_shape(graph: GeoGraph, source, limit,
                                       weight="length",
                                       alpha_quantile=99.0,
                                       remove_holes=True,
                                       tolerance=1e-7) -> GenericPolygon:
    """Returns an approximation of the isochrone polygon using an alpha-shape of the Shortest Path Tree.

    :param graph: GeoGraph to browse
    :param source: Source node from where distance is computed
    :param limit: Isochrone limit (e.g. 100 meters, 5 minutes, depending on ``weight`` unit).
    :param weight: Weight attribute on edges to compute distances (edge weights should be non-negative).
    :param alpha_quantile: Quantile on circumradius to determine alpha (100 returns the convex hull,
        0 returns an empty polygon). ``0 <= quantile <= 100``.
    :param remove_holes: If ``True`` remove holes in the returned polygon.
    :param tolerance: Buffering tolerance on polygon for rendering
    :return: A polygon approximating the isochrone.
    """
    # Compute the ego-graph
    ego_graph = gnx.extended_ego_graph(graph, source, limit, distance=weight)
    edge_as_lines = ego_graph.get_edges_as_line_series()
    discretized_lines, _ = gnx.discretize_lines(edge_as_lines)
    alpha_shape = get_alpha_shape_polygon(discretized_lines, alpha_quantile)
    alpha_shape = alpha_shape.buffer(tolerance)
    if remove_holes:
        if isinstance(alpha_shape, MultiPolygon):
            alpha_shape = MultiPolygon([Polygon(sub_p.exterior) for sub_p in alpha_shape])
        else:
            alpha_shape = Polygon(alpha_shape.exterior)
    return alpha_shape
