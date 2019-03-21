# -*- coding: utf-8 -*-
"""
    File name: test_simplify
    Author: Artelys
    Creation date: 23/01/2019
    Python Version: 3.6
"""

from nose.plugins.attrib import attr
from nose.tools import assert_less, assert_not_in, assert_in
import unittest
import geonetworkx as gnx
from geonetworkx.testing.utils import get_random_geograph_subclass, assert_is_subgraph, ALL_CLASSES
import geonetworkx.testing.utils as gnx_tu
from shapely.geometry import Polygon
import math
import numpy as np


NB_VERTICES = 50
gnx_tu.SEED = 70595

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

    def test_remove_nan_attributes(self):
        nan_examples = {'t1': None, 't2': math.nan, 't3': np.nan}
        not_nan_examples = {'t4': "abcd", 't5': 1.23, 't6': ["a", "b", "c", None, np.nan, math.nan]}
        for graph_type in ALL_CLASSES:
            for copy in [True, False]:
                with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED, copy=copy):
                    g = get_random_geograph_subclass(NB_VERTICES, graph_type)
                    for n, d in g.nodes(data=True):
                        d.update(nan_examples)
                        d.update(not_nan_examples)
                    for u, v, d in g.edges(data=True):
                        d.update(nan_examples)
                        d.update(not_nan_examples)
                    if copy:
                        modified_graph = gnx.remove_nan_attributes(g, copy=copy)
                    else:
                        gnx.remove_nan_attributes(g, copy=copy)
                        modified_graph = g

                    for n, d in modified_graph.nodes(data=True):
                        for k in nan_examples.keys():
                            assert_not_in(k, d.keys())
                        for k in not_nan_examples.keys():
                            assert_in(k, d.keys())
                    for u, v, d in modified_graph.edges(data=True):
                        for k in nan_examples.keys():
                            assert_not_in(k, d.keys())
                        for k in not_nan_examples.keys():
                            assert_in(k, d.keys())

