# -*- coding: utf-8 -*-
"""
    File name: isochrones
    Author: Artelys
    Creation date: 22/05/2019
    Python Version: 3.6
"""
import numpy as np
from geonetworkx.geograph import GeoGraph
import geonetworkx as gnx
from shapely.geometry import Polygon, LineString, MultiPolygon
from shapely.ops import cascaded_union
import math
from typing import Union
from scipy.spatial import Delaunay


def isochrone_subgraph(graph:GeoGraph, source, limit=500, weight_attr="length"):
    import networkx as nx
    dist, paths = nx.single_source_dijkstra(graph, source, cutoff=limit, weight=weight_attr)







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





def boundary_edge_buffer(line: LineString) -> Union[Polygon, MultiPolygon]:
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





def get_alpha_shape_polygon(points: list, quantile: int) -> Union[Polygon, MultiPolygon]:
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
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        # Semiperimeter of triangle
        s = (a + b + c) / 2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        circum_radius.append(circum_r)
        polygons.append(Polygon([pa, pb, pc]))
    alpha = np.percentile(circum_radius, quantile)
    filtered_polygons = [p for i, p in enumerate(polygons) if circum_radius[i] <= alpha]
    return cascaded_union(filtered_polygons)
