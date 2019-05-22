# -*- coding: utf-8 -*-
"""
    File name: isochrones
    Author: Artelys
    Creation date: 22/05/2019
    Python Version: 3.6
"""
import geonetworkx as gnx


# Read data set: Street graph sample of Grenoble, France.
graph = gnx.read_gpickle("../tests/datasets/grenoble_isere_streets.gpickle")

# Clean data: remove edges orientation, remove self-loops
graph = graph.to_undirected()
gnx.remove_self_loop_edges(graph)

# Supplement graph data: edges geometries, edges length
gnx.fill_edges_missing_geometry_attributes(graph)
gnx.fill_length_attribute(graph, attribute_name="length", only_missing=True)

# Setting isochrone parameters: source, limits
source = 2192254241
limits = [200, 400, 600]  # meters

# Compute isochrone polygon
polygons = []
for l in limits:
    p = gnx.isochrone_polygon(graph, source, l, weight="length")
    polygons.append(p)


