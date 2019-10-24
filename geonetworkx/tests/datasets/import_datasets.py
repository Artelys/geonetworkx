import os
import networkx as nx
import geonetworkx as gnx
import geopandas as gpd


def get_copenhagen_street_net():
    """Reads dataset geojson and return the corresponding geograph."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    nodes_path = os.path.join(dir_path, "copenhagen_streets_net_nodes.geojson")
    edges_path = os.path.join(dir_path, "copenhagen_streets_net_edges.geojson")
    return gnx.read_geofiles(nodes_path, edges_path, directed=True, multigraph=True)


def get_copenhagen_ferry_net():
    """Reads dataset geojson and return the corresponding geograph."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    nodes_path = os.path.join(dir_path, "copenhagen_ferry_net_nodes.geojson")
    edges_path = os.path.join(dir_path, "copenhagen_ferry_net_edges.geojson")
    return gnx.read_geofiles(nodes_path, edges_path, directed=True, multigraph=True)


def get_copenhagen_buildings():
    """Reads GeoDataFrame from geojson dataset"""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    return gpd.read_file(os.path.join(dir_path, "copenhagen_buildings.geojson"))


def get_grenoble_streets_500():
    """Reads dataset geojson and return the corresponding geograph."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    nodes_path = os.path.join(dir_path, "grenoble_streets_500_nodes.geojson")
    edges_path = os.path.join(dir_path, "grenoble_streets_500_edges.geojson")
    return gnx.read_geofiles(nodes_path, edges_path, directed=False, multigraph=True)


def get_grenoble_streets_200():
    """Reads dataset geojson and return the corresponding geograph."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    nodes_path = os.path.join(dir_path, "grenoble_streets_200_nodes.geojson")
    edges_path = os.path.join(dir_path, "grenoble_streets_200_edges.geojson")
    return gnx.read_geofiles(nodes_path, edges_path, directed=False, multigraph=True)


def get_grenoble_streets_isere():
    """Reads dataset geojson and return the corresponding geograph."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    nodes_path = os.path.join(dir_path, "grenoble_isere_nodes.geojson")
    edges_path = os.path.join(dir_path, "grenoble_isere_edges.geojson")
    return gnx.read_geofiles(nodes_path, edges_path, directed=True, multigraph=True)


def import_paris_subway_dataset():
    try:
        import osmnx as ox
    except ImportError:
        raise RuntimeError("OSMNX package is required for import data from OSM")
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
    def is_not_station(n):
        return gg.nodes[n].get("name", None) is None

    all_edges_names = nx.get_edge_attributes(gg, "name")
    gnx.remove_dead_ends(gg, node_filter=is_not_station, only_strict=False)
    edges_mapping = gnx.two_degree_node_merge(gg, node_filter=is_not_station)
    # Supplement data in the added edges with the most present value in the merged edges.
    for new_edge, removed_edges in edges_mapping.items():
        available_names = [all_edges_names[e] for e in removed_edges if e in all_edges_names]
        if available_names:
            gg.edges[new_edge]["name"] = max(set(available_names), key=available_names.count)

    # Save geograph
    dir_path = os.path.dirname(os.path.realpath(__file__))
    gnx.write_gpickle(gg, os.path.join(dir_path, "paris_subway.gpickle"))
