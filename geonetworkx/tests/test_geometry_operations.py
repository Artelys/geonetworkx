# -*- coding: utf-8 -*-
from shapely.geometry import LineString, MultiLineString, Point
import numpy as np
from nose.plugins.attrib import attr
import unittest
from geonetworkx.geometry_operations import discretize_lines
from geonetworkx.testing.utils import assert_almost_intersect
import geonetworkx as gnx
import geonetworkx.testing.utils as gnx_tu
import geonetworkx.settings as settings

gnx_tu.SEED = 70595
np.random.seed(gnx_tu.SEED)


@attr('geometry_operations')
class TestGeometryOperations(unittest.TestCase):

    def test_discretize_lines(self):
        # Test that a discretization of lines is well set for each line
        nb_lines = 20
        lines = []
        for l in range(nb_lines):
            line_nb_points = np.random.randint(5, 30)
            points = np.random.rand(line_nb_points, 2)
            lines.append(LineString(points))
        discretized_points, points_line_association = discretize_lines(lines, 1e-1)
        for l in points_line_association:
            line = lines[l]
            for p in points_line_association[l]:
                point = discretized_points[p]
                assert_almost_intersect(point, line, "A discretized point does not intersects its given line")

    def test_get_closest_point_from_multiline(self):
        l1 = LineString([(-73.615505, 45.504573),
                         (-73.614989, 45.504329),
                         (-73.614574, 45.504121),
                         (-73.614322, 45.503854),
                         (-73.614193, 45.503532),
                         (-73.614075, 45.503182),
                         (-73.613997, 45.503064),
                         (-73.613952, 45.503005)])
        l2 = LineString([(-73.615454, 45.504640),
                         (-73.614535, 45.504192),
                         (-73.614193, 45.503850),
                         (-73.613940, 45.503213),
                         (-73.613840, 45.503096)])
        ml = MultiLineString([l1, l2])
        points = [Point(-73.613941, 45.504743),
                  Point(-73.616263, 45.504509),
                  Point(-73.614920, 45.503820),
                  Point(-73.613956, 45.502349)]
        points_coords = np.array([[p.x, p.y] for p in points])
        d, ix = gnx.get_closest_point_from_multi_shape(ml, points_coords, discretization_tol=1e-5)
        self.assertEquals(ix, 2)
