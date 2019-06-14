# -*- coding: utf-8 -*-
"""
This script allows to download and formatting the Paris subway network to a geograph.
"""
import os
import networkx as nx
import osmnx as ox
import geonetworkx as gnx

# Set additional info to retrieve
ox.settings.useful_tags_node.append("name")
ox.settings.useful_tags_path.append("name")
# Download raw data with osmnx
g = ox.graph_from_address("Paris, France", distance=7000, infrastructure='relation["route"="subway"]',
                          clean_periphery=False, simplify=False, custom_filter="")

# Transform to a geograph
gg = gnx.read_geograph_with_coordinates_attributes(g)
gg.name = "Paris_subway"

# Supplement graph data
# Add geometry on all edges
gnx.fill_edges_missing_geometry_attributes(gg)
# Simplify: remove dead ends and two degree nodes that are not a subway station
is_not_station = lambda n: gg.nodes[n].get("name", None) is None
all_edges_names = nx.get_edge_attributes(gg, "name")
gnx.remove_dead_ends(gg, node_filter=is_not_station , only_strict=False)
edges_mapping = gnx.two_degree_node_merge(gg, node_filter=is_not_station)
# Supplement data in the added edges with the most present value in the merged edges.
for new_edge, removed_edges in edges_mapping.items():
    available_names = [all_edges_names[e] for e in removed_edges if e in all_edges_names]
    if available_names:
        gg.edges[new_edge]["name"] = max(set(available_names), key=available_names.count)


# Save geograph
dir_path = os.path.dirname(os.path.realpath(__file__))
gnx.write_gpickle(gg, os.path.join(dir_path, "paris_subway.gpickle"))