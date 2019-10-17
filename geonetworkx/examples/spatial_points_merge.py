# -*- coding: utf-8 -*-
""" This example shows how a two geographs can be merged in a routing use case. Here, a street network and a subway
network is combined to find the multi-modal shortest path between a source and a target."""
import os
import geopandas as gpd
import networkx as nx
import geonetworkx as gnx
import osmnx as ox


# Constants: speed in each network
bicycle_speed = 10 * 1000 / 3600  # 10 km/h -> 2.8 m/s

# Download and set up the street network (main streets_graph only)
streets_graph = ox.graph_from_address("Rennes, France", distance=2500,
                                      infrastructure='way["highway"~"primary|secondary|tertiary"]')
streets_graph = gnx.read_geograph_with_coordinates_attributes(streets_graph)
streets_graph.name = "Rennes_streets"

# Getting the bicycle stations
datasets_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../tests/datasets/")
bicycle_stations = gpd.read_file(os.path.join(datasets_path, "rennes_bicycle_stations_velo_star.geojson"))

# Merging the stations to the street network
gnx.spatial_points_merge(streets_graph, bicycle_stations, inplace=True)

# Setting times on the edges
gnx.fill_edges_missing_geometry_attributes(streets_graph)
gnx.fill_length_attribute(streets_graph, attribute_name="length", only_missing=True)
streets_lengths = nx.get_edge_attributes(streets_graph, "length")
nx.set_edge_attributes(streets_graph, {e: l / bicycle_speed for e, l in streets_lengths.items()}, "time")


# Computing all shortest path between all stations with a 15min cutoff
cutoff = 15 * 60  # in seconds
bicycle_stations_nodes = list(bicycle_stations.index)

# Creating a graph connecting all pairs of reachable stations within 15 min
bicycle_stations_graph = gnx.GeoDiGraph()
for n in bicycle_stations_nodes:
    bicycle_stations_graph.add_node(n, geometry=bicycle_stations.at[n, "geometry"])

for n in bicycle_stations_nodes:
    distances = nx.single_source_dijkstra_path_length(streets_graph,
                                                      source=n,
                                                      weight="time",
                                                      cutoff=cutoff)
    for u, d in distances.items():
        if u != n and u in bicycle_stations_nodes:
            bicycle_stations_graph.add_edge(n, u, time=d)

bicycle_stations_graph.name = "bicycle_station_graph"
gnx.write_geofile(bicycle_stations_graph, datasets_path)



