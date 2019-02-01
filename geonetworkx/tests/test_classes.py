"""
    File name: test_classes
    Author: Artelys
    Creation date: 17/01/2019
    Python Version: 3.6
"""
from geonetworkx.testing import get_random_geograph, get_random_geomultigraph, get_random_geodigraph,\
    get_random_geomultidigraph, get_random_geograph_with_wgs84_scale, get_random_geograph_subclass
from geonetworkx.testing import assert_graphs_have_same_edges_geometry, assert_graphs_have_same_geonodes, ALL_CLASSES
import geonetworkx as gnx
import os, shutil
from nose.tools import assert_is_instance, assert_equal
from nose.plugins.attrib import attr
import unittest

SEED = 70595
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
        graph = get_random_geograph(NB_POINTS)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoDiGraph)

    def test_graph_to_undirected(self):
        graph = get_random_geograph(NB_POINTS)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoGraph)

    # GeoMultiGraph
    def test_multigraph_to_directed(self):
        graph = get_random_geomultigraph(NB_POINTS)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoMultiDiGraph)

    def test_multigraph_to_undirected(self):
        graph = get_random_geomultigraph(NB_POINTS)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoMultiGraph)

    # GeoDiGraph
    def test_digraph_to_directed(self):
        graph = get_random_geodigraph(NB_POINTS)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoDiGraph)

    def test_digraph_to_undirected(self):
        graph = get_random_geodigraph(NB_POINTS)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoGraph)

    # GeoMultiDiGraph
    def test_multidigraph_to_directed(self):
        graph = get_random_geomultidigraph(NB_POINTS)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoMultiDiGraph)

    def test_multidigraph_to_undirected(self):
        graph = get_random_geomultidigraph(NB_POINTS)
        undirected_graph = graph.to_undirected()
        assert_is_instance(undirected_graph, gnx.GeoMultiGraph)

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
            g = get_random_geograph_subclass(NB_POINTS, graph_type)
            gdf = g.nodes_to_gdf()
            assert_equal(g.number_of_nodes(), len(gdf))

    def test_edges_to_gdf(self):
        for graph_type in ALL_CLASSES:
            g = get_random_geograph_subclass(NB_POINTS, graph_type)
            gdf = g.edges_to_gdf()
            assert_equal(g.number_of_edges(), len(gdf))

