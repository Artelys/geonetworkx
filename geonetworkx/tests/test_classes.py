"""
    File name: test_classes
    Author: Artelys
    Creation date: 17/01/2019
    Python Version: 3.6
"""
from geonetworkx.testing import get_random_geograph, get_random_geomultigraph,\
                                get_random_geodigraph, get_random_geomultidigraph
import geonetworkx as gnx
import numpy as np
from nose.tools import assert_is_instance


SEED = 70595
np.random.seed(SEED)
NB_POINTS = 100

class TestClasses():

    # GeoGraph
    def test_graph_to_directed(self):
        graph = get_random_geograph(NB_POINTS, SEED)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoDiGraph)

    def test_graph_to_undirected(self):
        graph = get_random_geograph(NB_POINTS, SEED + 1)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoGraph)

    # GeoMultiGraph
    def test_multigraph_to_directed(self):
        graph = get_random_geomultigraph(NB_POINTS, SEED + 2)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoMultiDiGraph)

    def test_multigraph_to_undirected(self):
        graph = get_random_geomultigraph(NB_POINTS, SEED + 3)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoMultiGraph)

    # GeoDiGraph
    def test_digraph_to_directed(self):
        graph = get_random_geodigraph(NB_POINTS, SEED + 4)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoDiGraph)

    def test_digraph_to_undirected(self):
        graph = get_random_geodigraph(NB_POINTS, SEED + 5)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoGraph)

    # GeoMultiDiGraph
    def test_multidigraph_to_directed(self):
        graph = get_random_geomultidigraph(NB_POINTS, SEED + 6)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoMultiDiGraph)

    def test_multidigraph_to_undirected(self):
        graph = get_random_geomultidigraph(NB_POINTS, SEED + 7)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoMultiGraph)




