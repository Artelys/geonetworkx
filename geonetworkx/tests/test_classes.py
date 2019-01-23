"""
    File name: test_classes
    Author: Artelys
    Creation date: 17/01/2019
    Python Version: 3.6
"""
from geonetworkx.testing import get_random_geograph, get_random_geomultigraph, get_random_geodigraph,\
    get_random_geomultidigraph, get_random_geograph_with_wgs84_scale
from geonetworkx.testing import assert_graphs_have_same_edges_geometry, assert_graphs_have_same_geonodes
import geonetworkx as gnx
import os, shutil
import numpy as np
from nose.tools import assert_is_instance
from nose.plugins.attrib import attr
import unittest


SEED = 70595
np.random.seed(SEED)
NB_POINTS = 100

@attr('classes')
class TestClasses(unittest.TestCase):

    def setUp(self):
        file_dir = os.path.dirname(__file__)
        self.results_dir = os.path.join(file_dir, "datasets/results/")
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def tearDown(self):
        shutil.rmtree(self.results_dir)

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

    def test_crs_modification(self):
        graph = get_random_geograph_with_wgs84_scale(NB_POINTS, SEED, gnx.GeoMultiDiGraph)
        modified_graph = graph.to_crs(crs={'init': 'epsg:3857'}, inplace=False)
        re_modified_graph = modified_graph.to_crs(crs=gnx.settings.WGS84_CRS, inplace=False)
        assert_graphs_have_same_edges_geometry(graph, re_modified_graph, "Some edge geometries seems to be different"
                                                                         " after re-modification", tol=1e-2)
        assert_graphs_have_same_geonodes(graph, re_modified_graph, "Some nodes seems to be different after "
                                                                   "re-modification")





