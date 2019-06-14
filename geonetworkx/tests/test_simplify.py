# -*- coding: utf-8 -*-
"""
    File name: test_simplify
    Author: Artelys
    Creation date: 23/01/2019
    Python Version: 3.6
"""

from nose.plugins.attrib import attr
from nose.tools import assert_less, assert_not_in, assert_in, assert_equal, assert_true, assert_false
import unittest
import geonetworkx as gnx
from geonetworkx.testing.utils import get_random_geograph_subclass, assert_is_subgraph, ALL_CLASSES
import geonetworkx.testing.utils as gnx_tu
from shapely.geometry import Polygon, Point
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

    def test_two_degree_merge(self):
        nodes_coordinates = [(1, {"geometry": Point(0, 1)}), (2, {"geometry": Point(1, 1)}),
                             (3, {"geometry": Point(2, 1)}), (4, {"geometry": Point(1, 0)}),
                             (5, {"geometry": Point(1, 2)}), (6, {"geometry": Point(3, 1)}),
                             (7, {"geometry": Point(1, 3)}), (8, {"geometry": Point(1, 4)}),
                             (9, {"geometry": Point(1, 5)}), (10, {"geometry": Point(2, 5)}),
                             (11, {"geometry": Point(3, 5)}), (12, {"geometry": Point(2, 6)}),
                             (13, {"geometry": Point(3, 6)})]
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type):
                g = graph_type()
                g.add_nodes_from(nodes_coordinates)
                g.add_path([1, 2, 3])
                g.add_path([4, 2, 5])
                g.add_path([3, 6])
                g.add_path([5, 7, 8, 9])
                g.add_path([11, 10, 9])
                g.add_path([9, 12, 13])
                if g.is_directed():
                    g.add_path([6, 3])
                    g.add_path([13, 12, 9])
                    if g.is_multigraph():
                        g.add_edge(10, 9)
                        merged_nodes = {5, 7, 8, 12}
                        added_edges = {(2, 9), (9, 13), (13, 9)}
                    else:
                        merged_nodes = {5, 7, 8, 10, 12}
                        added_edges = {(2, 9), (11, 9), (9, 13), (13, 9)}
                else:
                    if g.is_multigraph():
                        g.add_edge(10, 9)
                        merged_nodes = {3, 5, 7, 8, 12}
                        added_edges = {(2, 9), (2, 6)}
                    else:
                        merged_nodes = {3, 5, 7, 8, 10, 12}
                        added_edges = {(2, 9), (2, 6), (9, 11), (9, 13)}
                initial_nb_nodes = g.number_of_nodes()
                gnx.fill_edges_missing_geometry_attributes(g)
                if g.is_multigraph():
                    e = (5, 7, 0)
                else:
                    e = (5, 7)
                del g.edges[e][g.edges_geometry_key]
                old_g = g.copy()
                merged_edges = gnx.two_degree_node_merge(g)
                for n in merged_nodes:
                    assert_not_in(n, g.nodes(), "A merged node is a in the simplified graph: %s" % str(n))
                assert_equal(g.number_of_nodes(), initial_nb_nodes - len(merged_nodes))
                for e in added_edges:
                    assert_true(g.has_edge(*e), "A normally added edge is not in the simplified graph: %s" % str(e))
                for new_edge, edges in merged_edges.items():
                    assert_true(g.has_edge(*new_edge), "A new edge is not in the simplified graph: %s" % str(e))
                    for e in edges:
                        assert_false(g.has_edge(*e), "A deleted edge is in the simplified graph: %s" % str(e))
                        assert_true(old_g.has_edge(*e), "A deleted edge is not not in the old graph: %s" % str(e))

