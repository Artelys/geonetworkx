import os
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