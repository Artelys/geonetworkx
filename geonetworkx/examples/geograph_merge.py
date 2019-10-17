# -*- coding: utf-8 -*-
""" This example shows how a two geographs can be merged in a routing use case. Here, a street network and a subway
network is combined to find the multi-modal shortest path between a source and a target."""
import os
import networkx as nx
import geonetworkx as gnx
import osmnx as ox
from shapely.geometry import box

# Constants: speed in each network
bike_speed = 10 * 1000 / 3600  # 10 km/h -> 2.8 m/s
subway_speed = 50 * 1000 / 3600  # 50 km/h -> 13.9 m/s

# Download and set up the street network (this can be long depending on the connection speed)
streets = ox.graph_from_address("Paris, France", distance=2000)
streets = gnx.read_geograph_with_coordinates_attributes(streets)
streets.name = "Paris_streets"
# Setting times on the edges
streets_lengths = nx.get_edge_attributes(streets, "length")
nx.set_edge_attributes(streets, {e: l / bike_speed for e, l in streets_lengths.items()}, "time")

# Getting the subway network
dir_path = os.path.dirname(os.path.realpath(__file__))
subway = gnx.read_gpickle(os.path.join(dir_path, "paris_subway.gpickle"))
# Setting the times on the edges
gnx.fill_length_attribute(subway)
subway_lengths = nx.get_edge_attributes(subway, "length")
nx.set_edge_attributes(subway, {e: l / subway_speed for e, l in subway_lengths.items()}, "time")
bb = gnx.get_graph_bounding_box(streets)
streets_bb = box(bb[0][0], bb[0][1], bb[1][0], bb[1][1])


def subway_node_is_reachable_from_streets(n):
    """A subway node is reachable from streets if it is a station (i.e. it has a name) and if it is in the bounding box
    of the the street network."""
    if subway.nodes[n].get("name", None) is None:
        return False
    if not streets_bb.contains(subway.nodes[n][subway.nodes_geometry_key]):
        return False
    return True


# Merging the two graphs
streets_and_subway = gnx.spatial_graph_merge(streets, subway, node_filter=subway_node_is_reachable_from_streets)
streets_and_subway.name = "Paris_streets_subway"

# Computing the shortest path
path = nx.shortest_path(streets_and_subway, source=25273868, target=4518793163)
# Saving the optimal path on edges with an attribute
edge_like_path = [(path[i], path[i + 1], 0) for i in range(len(path) - 1)]
for e in edge_like_path:
    streets_and_subway.edges[e]["opt_path"] = 1



