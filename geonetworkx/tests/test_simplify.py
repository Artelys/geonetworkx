"""
    File name: test_simplify
    Author: Artelys
    Creation date: 23/01/2019
    Python Version: 3.6
"""

from nose.plugins.attrib import attr
from nose.tools import assert_less
import unittest
import geonetworkx as gnx
from geonetworkx.testing.utils import get_random_geograph_subclass, assert_is_subgraph, ALL_CLASSES
from shapely.geometry import Polygon


NB_VERTICES = 50
SEED = 70595

@attr("simplify")
class TestSimplify(unittest.TestCase):

    def test_trim_graph_with_polygon(self):
        for method in ["intersects", "within"]:
            for graph_type in ALL_CLASSES:
                g = get_random_geograph_subclass(NB_VERTICES, graph_type)
                nodes_series = g.get_nodes_as_point_series()
                bounds = nodes_series.total_bounds
                bounds[1] = (bounds[1] + bounds[3]) / 2.0  # Taking half bounding box as polygon
                polygon_coordinates = [[bounds[0], bounds[1]],
                                        [bounds[0], bounds[3]],
                                        [bounds[2], bounds[3]],
                                        [bounds[2], bounds[1]]]
                polygon = Polygon(polygon_coordinates)
                new_graph = gnx.trim_graph_with_polygon(g, polygon, copy=True, method=method)
                assert_less(len(new_graph.nodes), len(g.nodes), "Half bounding box trimming must remove some nodes")
                assert_is_subgraph(g, new_graph, "The trimmed graph must be a sub graph of the original graph")
