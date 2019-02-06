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
import geonetworkx as gnx
from geonetworkx.testing import assert_graphs_have_same_edges_geometry, ALL_CLASSES, get_random_geograph_subclass
import geonetworkx.testing.utils as gnx_tu


gnx_tu.SEED = 70595
NB_VERTICES = 50


@attr('readwrite')
class TestReadWrite(unittest.TestCase):

    def setUp(self):
        file_dir = os.path.dirname(__file__)
        self.results_dir = os.path.join(file_dir, "datasets/results/")
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def tearDown(self):
        shutil.rmtree(self.results_dir)

    def test_write_read_gpickle(self):
        for i, graph_type in enumerate(ALL_CLASSES):
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_VERTICES, graph_type)
                file_path = os.path.join(self.results_dir, "write_test_%d.gpickle" % i)
                gnx.write_gpickle(g, file_path)
                read_graph = gnx.read_gpickle(file_path)
                assert_is_instance(read_graph, graph_type)
                assert_graphs_have_same_edges_geometry(g, read_graph)

    def test_write_read_graphml(self):
        for i, graph_type in enumerate(ALL_CLASSES):
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_VERTICES, graph_type)
                gnx.stringify_nodes(g, copy=False)
                file_path = os.path.join(self.results_dir, "write_test_%d.xml" % i)
                gnx.write_graphml(g, file_path)
                read_graph = gnx.read_graphml(file_path)
                assert_is_instance(read_graph, graph_type)
                assert_graphs_have_same_edges_geometry(g, read_graph)


    def test_write_geofile(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_VERTICES, graph_type)
                gnx.write_geofile(g, self.results_dir, driver="GPKG")

