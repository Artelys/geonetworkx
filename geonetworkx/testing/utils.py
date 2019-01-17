"""
    File name: utils
    Author: Artelys
    Creation date: 08/01/2019
    Python Version: 3.6
"""
import networkx as nx
import geonetworkx as gnx
import numpy as np
from nose.tools import assert_true


def assert_almost_intersect(shape1, shape2, msg='', tol=1e-4):
    assert_true(shape1.buffer(tol).intersects(shape2), msg)


def get_random_geograph(nb_nodes, seed):
    edge_creation_prob = 0.1
    g = nx.fast_gnp_random_graph(nb_nodes, edge_creation_prob, seed=seed, directed=False)
    nodes_coords = nx.kamada_kawai_layout(g)
    nx.set_node_attributes(g, {n: coords[0] for n, coords in nodes_coords.items()}, 'x')
    nx.set_node_attributes(g, {n: coords[1] for n, coords in nodes_coords.items()}, 'y')
    graph = gnx.GeoGraph(g)
    gnx.fill_edges_missing_geometry_attributes(graph)
    return graph

def get_random_geodigraph(nb_nodes, seed):
    edge_creation_prob = 0.1
    g = nx.fast_gnp_random_graph(nb_nodes, edge_creation_prob, seed=seed, directed=True)
    nodes_coords = nx.kamada_kawai_layout(g)
    nx.set_node_attributes(g, {n: coords[0] for n, coords in nodes_coords.items()}, 'x')
    nx.set_node_attributes(g, {n: coords[1] for n, coords in nodes_coords.items()}, 'y')
    graph = gnx.GeoDiGraph(g)
    gnx.fill_edges_missing_geometry_attributes(graph)
    return graph

def get_random_geomultigraph(nb_nodes, seed):
    edge_creation_prob = 0.1
    g = get_random_geograph(nb_nodes, seed)
    multigraph = gnx.GeoMultiGraph(g)
    np.random.seed(seed)
    for e in list(multigraph.edges):
        if np.random.rand() < edge_creation_prob:
            multigraph.add_edge(e[0], e[1])
    gnx.fill_edges_missing_geometry_attributes(multigraph)
    return multigraph

def get_random_geomultidigraph(nb_nodes, seed):
    edge_creation_prob = 0.1
    g = get_random_geodigraph(nb_nodes, seed)
    multidigraph = gnx.GeoMultiGraph(g)
    np.random.seed(seed)
    for e in list(multidigraph.edges):
        if np.random.rand() < edge_creation_prob:
            multidigraph.add_edge(e[0], e[1])
    gnx.fill_edges_missing_geometry_attributes(multidigraph)
    return multidigraph