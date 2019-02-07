"""
    File name: test_classes
    Author: Artelys
    Creation date: 17/01/2019
    Python Version: 3.6
"""
from geonetworkx.testing import get_random_geograph_with_wgs84_scale, get_random_geograph_subclass
from geonetworkx.testing import assert_graphs_have_same_edges_geometry, assert_graphs_have_same_geonodes, ALL_CLASSES, \
                                assert_graphs_have_same_spatial_keys
import geonetworkx.testing.utils as gnx_tu
import geonetworkx as gnx
import os, shutil
from nose.tools import assert_is_instance, assert_equal
from nose.plugins.attrib import attr
import unittest

gnx_tu.SEED = 70595
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
        graphs_directed_match = {gnx.GeoGraph: gnx.GeoDiGraph, gnx.GeoMultiGraph: gnx.GeoMultiDiGraph,
                                 gnx.GeoDiGraph: gnx.GeoDiGraph, gnx.GeoMultiDiGraph: gnx.GeoMultiDiGraph}
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                graph = get_random_geograph_subclass(NB_POINTS, graph_type)
                directed_graph = graph.to_directed()
                assert_is_instance(directed_graph, graphs_directed_match[graph_type])

    def test_graph_to_undirected(self):
        graphs_undirected_match = {gnx.GeoGraph: gnx.GeoGraph, gnx.GeoMultiGraph: gnx.GeoMultiGraph,
                                   gnx.GeoDiGraph: gnx.GeoGraph, gnx.GeoMultiDiGraph: gnx.GeoMultiGraph}
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                graph = get_random_geograph_subclass(NB_POINTS, graph_type)
                directed_graph = graph.to_undirected()
                assert_is_instance(directed_graph, graphs_undirected_match[graph_type])

    def test_crs_modification(self):
        graph = get_random_geograph_with_wgs84_scale(NB_POINTS, gnx.GeoMultiDiGraph)
        modified_graph = graph.to_crs(crs={'init': 'epsg:3857'}, inplace=False)
        re_modified_graph = modified_graph.to_crs(crs=gnx.settings.WGS84_CRS, inplace=False)
        assert_graphs_have_same_edges_geometry(graph, re_modified_graph, "Some edge geometries seems to be different"
                                                                         " after re-modification", tol=1e-2)
        assert_graphs_have_same_geonodes(graph, re_modified_graph, "Some nodes seems to be different after "
                                                                   "re-modification")

    def test_nodes_to_gdf(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_POINTS, graph_type)
                gdf = g.nodes_to_gdf()
                assert_equal(g.number_of_nodes(), len(gdf))

    def test_edges_to_gdf(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_POINTS, graph_type)
                gdf = g.edges_to_gdf()
                assert_equal(g.number_of_edges(), len(gdf))

    def test_spatial_keys_persistence(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_POINTS, graph_type)
                gnx.rename_edges_attribute(g, g.edges_geometry_key, "abcd")
                g.edges_geometry_key = "abcd"
                gnx.rename_nodes_attribute(g, g.x_key, "efgh")
                g.x_key = "efgh"
                gnx.rename_nodes_attribute(g, g.y_key, "ijkl")
                g.y_key = "ijkl"
                g.crs = {'init': 'epsg:3945'}
                g2 = g.copy(as_view=False)
                assert_graphs_have_same_spatial_keys(g, g2)
                g3 = g.to_undirected(as_view=False)
                assert_graphs_have_same_spatial_keys(g, g3)
                g4 = g.to_directed(as_view=False)
                assert_graphs_have_same_spatial_keys(g, g4)



