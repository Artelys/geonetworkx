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


gnx_tu.SEED = 70595
NB_POINTS = 50

@attr('tools')
class TestUtils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_edges_lines_splitting(self):
        g = gnx_tu.get_random_geograph_with_wgs84_scale(NB_POINTS, gnx.GeoMultiDiGraph)
        lines = list(g.get_edges_as_line_series())
        import sys, os
        from shapely.geometry import *
        os.chdir(r"D:\projets\GeoNetworkX\GeoNetworkX\geonetworkx\utils")
        sys.path.append(os.getcwd())
        from voronoi_parser import PyVoronoiHelper
        res = PyVoronoiHelper.splitting_test(lines)

