"""
    File name: test_geometry_operations.py
    Author: Artelys
    Creation date: 08/01/2019
    Python Version: 3.6
"""
from shapely.geometry import LineString
import numpy as np
from nose.plugins.attrib import attr
import unittest
from geonetworkx.geometry_operations import discretize_lines
from geonetworkx.testing.utils import assert_almost_intersect
import geonetworkx.testing.utils as gnx_tu
import geonetworkx as gnx
import geonetworkx.settings as settings

gnx_tu.SEED = 70595
np.random.seed(gnx_tu.SEED)

@attr('geometry_operations')
class TestGeometryOperations(unittest.TestCase):

    def test_discretize_lines(self):
        # Test that a discretization of lines is well set for each line
        settings.DISCRETIZATION_TOLERANCE = 1e-1
        nb_lines = 20
        lines = []
        for l in range(nb_lines):
            line_nb_points = np.random.randint(5, 30)
            points = np.random.rand(line_nb_points, 2)
            lines.append(LineString(points))
        discretized_points, points_line_association = discretize_lines(lines)
        for l in points_line_association:
            line = lines[l]
            for p in points_line_association[l]:
                point = discretized_points[p]
                assert_almost_intersect(point, line, "A discretized point does not intersects its given line")


