"""
    File name: utils
    Author: Artelys
    Creation date: 08/01/2019
    Python Version: 3.6
"""
import networkx as nx
import geonetworkx as gnx
import numpy as np
from nose.tools import assert_true, assert_equal, assert_in


def assert_almost_intersect(shape1, shape2, msg='', tol=1e-4):
    assert_true(shape1.buffer(tol).intersects(shape2), msg)

def assert_coordinates_almost_equals(coordinates1, coordinates2, msg='', tol=1e-4):
    assert_true(gnx.coordinates_almost_equal(coordinates1, coordinates2, tol), msg)

def assert_lines_almost_equals(line1, line2, msg='', tol=1e-4):
    assert_equal(len(line1.coords), len(line2.coords), msg)
    for c1, c2 in zip(line1.coords, line2.coords):
        assert_coordinates_almost_equals(c1, c2, msg, tol)

def assert_graphs_have_same_edges_geometry(graph1, graph2, msg='', tol=1e-4):
    assert_equal(len(graph1.edges), len(graph2.edges), msg)
    lines1 = nx.get_edge_attributes(graph1, graph1.edges_geometry_key)
    lines2 = nx.get_edge_attributes(graph2, graph2.edges_geometry_key)
    for e in lines1:
        assert_in(e, lines2, msg)
        assert_lines_almost_equals(lines1[e], lines2[e], msg, tol)




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
    multidigraph = gnx.GeoMultiDiGraph(g)
    np.random.seed(seed)
    for e in list(multidigraph.edges):
        if np.random.rand() < edge_creation_prob:
            multidigraph.add_edge(e[0], e[1])
    gnx.fill_edges_missing_geometry_attributes(multidigraph)
    return multidigraph