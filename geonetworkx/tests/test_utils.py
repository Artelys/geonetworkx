# -*- coding: utf-8 -*-
import sys
import networkx as nx
import geonetworkx as gnx
import unittest
from nose.tools import assert_true, assert_equal, assert_less_equal
from nose.plugins.attrib import attr
import geonetworkx.testing.utils as gnx_tu
from geonetworkx.utils.voronoi_utils import *
from geonetworkx.generators import extended_ego_graph
from geonetworkx.tests import datasets
import os
import warnings


gnx_tu.SEED = 70595
NB_POINTS = 50
data_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), "datasets")


@attr('tools')
class TestUtils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_voronoi_edges(self):
        if not gnx_tu.check_optional_package_presence("pyvoronoi"):
            warnings.warn("Missing optional package for test 'test_voronoi_edges': pyvoronoi")
            return
        if sys.version_info.major >= 3 and sys.version_info.minor >= 8:
            warnings.warn("'pyvoronoi' package is not compatible with python 3.8 or later.")
            return
        gmg = datasets.get_grenoble_streets_200()
        gnx.fill_edges_missing_geometry_attributes(gmg)
        edge_as_lines = gmg.get_edges_as_line_series()
        lines = list(edge_as_lines)
        tolerance = 1e-7
        res = compute_voronoi_cells_from_lines(lines, tolerance)
        for e, line in edge_as_lines.items():
            cell_found = False
            for p in res:
                if p.buffer(10 * tolerance).contains(line):
                    cell_found = True
                    break
            if not cell_found:
                assert_true(False, "A edge geometry '%s' is not in any voronoi cells" % str(e))

    def test_extented_ego_graph(self):
        gmg = datasets.get_grenoble_streets_200()
        gnx.fill_edges_missing_geometry_attributes(gmg)
        gnx.fill_length_attribute(gmg, "length", only_missing=True)
        source = 312173744
        limit = 100  # meters
        ego_graph = extended_ego_graph(gmg, source, radius=limit, center=True, undirected=False, distance="length")
        ccs = list(nx.connected_components(ego_graph))
        assert_equal(len(ccs), 1, "The ego graph has several connected components")
        path_lengths = nx.single_source_dijkstra_path_length(ego_graph, source, weight="length")
        for n, length in path_lengths.items():
            assert_less_equal(length, limit, "A path in the ego graph is too long: %f > %f" % (length, limit))

