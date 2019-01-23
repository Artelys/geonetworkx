"""
    File name: test_readwrite.py
    Author: Artelys
    Creation date: 21/01/2019
    Python Version: 3.6
"""
import os, shutil
from nose.tools import assert_is_instance
import unittest
from nose.plugins.attrib import attr
import networkx as nx
import geonetworkx as gnx
from geonetworkx.testing import get_random_geograph, get_random_geodigraph, get_random_geomultigraph,\
    get_random_geomultidigraph, assert_graphs_have_same_edges_geometry


SEED = 70595
NB_VERTICES = 100

@attr('readwrite')
class TestReadWrite(unittest.TestCase):

    def setUp(self):
        file_dir = os.path.dirname(__file__)
        self.results_dir = os.path.join(file_dir, "datasets/results/")
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def tearDown(self):
        shutil.rmtree(self.results_dir)

    def test_write_read_gpickle_gg(self):
        g = get_random_geograph(NB_VERTICES, SEED)
        file_path = os.path.join(self.results_dir, "write_test_gg.gpickle")
        gnx.write_gpickle(g, file_path)
        read_graph = gnx.read_gpickle(file_path)
        assert_is_instance(read_graph, gnx.GeoGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)

    def test_write_read_gpickle_gdg(self):
        g = get_random_geodigraph(NB_VERTICES, SEED)
        file_path = os.path.join(self.results_dir, "write_test_gdg.gpickle")
        gnx.write_gpickle(g, file_path)
        read_graph = gnx.read_gpickle(file_path)
        assert_is_instance(read_graph, gnx.GeoDiGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)

    def test_write_read_gpickle_gmg(self):
        g = get_random_geomultigraph(NB_VERTICES, SEED)
        file_path = os.path.join(self.results_dir, "write_test_gmg.gpickle")
        gnx.write_gpickle(g, file_path)
        read_graph = gnx.read_gpickle(file_path)
        assert_is_instance(read_graph, gnx.GeoMultiGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)

    def test_write_read_gpickle_gmdg(self):
        g = get_random_geomultidigraph(NB_VERTICES, SEED)
        file_path = os.path.join(self.results_dir, "write_test_gmdg.gpickle")
        gnx.write_gpickle(g, file_path)
        read_graph = gnx.read_gpickle(file_path)
        assert_is_instance(read_graph, gnx.GeoMultiDiGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)


    def test_write_read_graphml_gg(self):
        g = get_random_geograph(NB_VERTICES, SEED)
        gnx.stringify_nodes(g, copy=False)
        file_path = os.path.join(self.results_dir, "write_test_gg.xml")
        gnx.write_graphml(g, file_path)
        read_graph = gnx.read_graphml(file_path)
        assert_is_instance(read_graph, gnx.GeoGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)

    def test_write_read_graphml_gdg(self):
        g = get_random_geodigraph(NB_VERTICES, SEED)
        gnx.stringify_nodes(g, copy=False)
        file_path = os.path.join(self.results_dir, "write_test_gdg.xml")
        gnx.write_graphml(g, file_path)
        read_graph = gnx.read_graphml(file_path)
        assert_is_instance(read_graph, gnx.GeoDiGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)

    def test_write_read_graphml_gmg(self):
        g = get_random_geomultigraph(NB_VERTICES, SEED)
        gnx.stringify_nodes(g, copy=False)
        file_path = os.path.join(self.results_dir, "write_test_gmg.xml")
        gnx.write_graphml(g, file_path)
        read_graph = gnx.read_graphml(file_path)
        assert_is_instance(read_graph, gnx.GeoMultiGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)

    def test_write_read_graphml_gmdg(self):
        g = get_random_geomultidigraph(NB_VERTICES, SEED)
        gnx.stringify_nodes(g, copy=False)
        file_path = os.path.join(self.results_dir, "write_test_gmdg.xml")
        gnx.write_graphml(g, file_path)
        read_graph = gnx.read_graphml(file_path)
        assert_is_instance(read_graph, gnx.GeoMultiDiGraph)
        assert_graphs_have_same_edges_geometry(g, read_graph)

