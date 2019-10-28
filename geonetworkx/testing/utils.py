# -*- coding: utf-8 -*-
import networkx as nx
import geonetworkx as gnx
import numpy as np
from shapely.geometry import LineString, Point
from nose.tools import assert_true, assert_equal, assert_in
from geonetworkx.utils import crs_equals
import importlib


ALL_CLASSES = [gnx.GeoGraph, gnx.GeoMultiGraph, gnx.GeoDiGraph, gnx.GeoMultiDiGraph]
SEED = 70595


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


def assert_is_subgraph(base_graph, sub_graph, msg=''):
    for n in sub_graph.nodes:
        assert_in(n, base_graph.nodes, msg)
    for e in sub_graph.edges:
        assert_in(e, base_graph.edges, msg)


def assert_graphs_have_same_geonodes(graph1, graph2, msg='', tol=1e-4):
    assert_equal(len(graph1.nodes), len(graph2.nodes), msg)
    coordinates1 = graph1.get_nodes_coordinates()
    coordinates2 = graph1.get_nodes_coordinates()
    assert_equal(len(coordinates1), len(coordinates2), msg)
    for n in coordinates1.keys():
        assert_in(n, coordinates2)
        assert_coordinates_almost_equals(coordinates1[n], coordinates2[n], msg, tol)


def assert_graphs_have_same_spatial_keys(graph1, graph2, msg=''):
    keys1 = graph1.get_spatial_keys()
    keys2 = graph2.get_spatial_keys()
    assert_equal(len(keys1), len(keys2))
    for k in keys1.keys():
        if k != 'crs':
            assert_equal(keys1[k], keys2[k], msg)
        else:
            assert_true(crs_equals(keys1[k], keys2[k]), msg)


def get_random_geograph(nb_nodes):
    global SEED
    edge_creation_prob = 0.1
    g = nx.fast_gnp_random_graph(nb_nodes, edge_creation_prob, seed=SEED, directed=False)
    SEED += 1
    nodes_coords = nx.kamada_kawai_layout(g)
    geometry_attr = "geometry"
    nx.set_node_attributes(g, {n: Point(coords) for n, coords in nodes_coords.items()}, geometry_attr)
    graph = gnx.GeoGraph(g, nodes_geometry_key=geometry_attr, edges_geometry_key=geometry_attr)
    gnx.fill_edges_missing_geometry_attributes(graph)
    return graph


def get_random_geodigraph(nb_nodes):
    global SEED
    edge_creation_prob = 0.1
    g = nx.fast_gnp_random_graph(nb_nodes, edge_creation_prob, seed=SEED, directed=True)
    SEED += 1
    nodes_coords = nx.kamada_kawai_layout(g)
    geometry_attr = "geometry"
    nx.set_node_attributes(g, {n: Point(coords) for n, coords in nodes_coords.items()}, geometry_attr)
    graph = gnx.GeoDiGraph(g, nodes_geometry_key=geometry_attr, edges_geometry_key=geometry_attr)
    gnx.fill_edges_missing_geometry_attributes(graph)
    return graph


def get_random_geomultigraph(nb_nodes):
    global SEED
    edge_creation_prob = 0.1
    g = get_random_geograph(nb_nodes)
    multigraph = gnx.GeoMultiGraph(g)
    np.random.seed(SEED)
    SEED += 1
    for e in list(multigraph.edges):
        if np.random.rand() < edge_creation_prob:
            multigraph.add_edge(e[0], e[1])
    gnx.fill_edges_missing_geometry_attributes(multigraph)
    return multigraph


def get_random_geomultidigraph(nb_nodes):
    global SEED
    edge_creation_prob = 0.1
    g = get_random_geodigraph(nb_nodes)
    multidigraph = gnx.GeoMultiDiGraph(g)
    np.random.seed(SEED)
    SEED += 1
    for e in list(multidigraph.edges):
        if np.random.rand() < edge_creation_prob:
            multidigraph.add_edge(e[0], e[1])
    gnx.fill_edges_missing_geometry_attributes(multidigraph)
    return multidigraph


def get_random_geograph_subclass(nb_nodes, graph_type=gnx.GeoGraph):
    graph_generation_methods = {gnx.GeoGraph: get_random_geograph,
                                gnx.GeoMultiGraph: get_random_geomultigraph,
                                gnx.GeoDiGraph: get_random_geodigraph,
                                gnx.GeoMultiDiGraph: get_random_geomultidigraph}
    g = graph_generation_methods[graph_type](nb_nodes)
    return g


def get_random_geograph_with_wgs84_scale(nb_nodes, graph_type=gnx.GeoGraph):
    g = get_random_geograph_subclass(nb_nodes, graph_type)
    x_func = lambda x: 5.7 + 1e-1 * x
    y_func = lambda y: 45.1 + 1e-1 * y
    transform_point = lambda p: Point([x_func(p.x), y_func(p.y)])
    for n, data in g.nodes(data=True):
        data[g.nodes_geometry_key] = transform_point(data[g.nodes_geometry_key])
    for e in g.edges:
        edge_data = g.edges[e]
        line = edge_data[g.edges_geometry_key]
        line_coords = list(line.coords)
        modified_coords = [(x_func(x), y_func(y)) for x, y in line_coords]
        edge_data[g.edges_geometry_key] = LineString(modified_coords)
    g.crs = gnx.settings.WGS84_CRS
    return g


def check_optional_package_presence(pkg_name):
    """Returns true if the given package name is available, false otherwise."""
    try:
        importlib.import_module(pkg_name)
    except ImportError as e:
        if pkg_name in e.msg:
            return False
    return True
