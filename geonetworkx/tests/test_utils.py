# -*- coding: utf-8 -*-
"""
    File name: test_utils
    Author: Artelys
    Creation date: 06/05/2019
    Python Version: 3.6
"""

import networkx as nx
import geonetworkx as gnx
import unittest
from nose.plugins.attrib import attr
import geonetworkx.testing.utils as gnx_tu
from geonetworkx.utils import voronoi_parser
import os


gnx_tu.SEED = 70595
NB_POINTS = 50
data_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), "datasets")


@attr('tools')
class TestUtils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_edges_lines_splitting(self):
        g = gnx_tu.get_random_geograph_with_wgs84_scale(NB_POINTS, gnx.GeoMultiDiGraph)
        lines = list(g.get_edges_as_line_series())
        res = voronoi_parser.compute_voronoi_cells_from_lines(lines)

    def test_voronoi_edges(self):
        mdg = nx.read_gpickle(os.path.join(data_directory, "grenoble200_mdg.gpickle"))
        mg = mdg.to_undirected()
        gmg = gnx.read_geograph_with_coordinates_attributes(mg)
        lines = list(gmg.get_edges_as_line_series())
        res = compute_voronoi_cells_from_lines(lines)
        res.to_file(r"C:\Users\hchareyre\Documents\trash\Nouveau dossier (8)\test_gre.shp")


