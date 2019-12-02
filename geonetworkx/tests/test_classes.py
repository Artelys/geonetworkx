# -*- coding: utf-8 -*-
from geonetworkx.testing import get_random_geograph_with_wgs84_scale, get_random_geograph_subclass
from geonetworkx.testing import assert_graphs_have_same_edges_geometry, assert_graphs_have_same_geonodes, ALL_CLASSES, \
                                assert_graphs_have_same_spatial_keys, assert_points_almost_equals
import geonetworkx.testing.utils as gnx_tu
import geonetworkx as gnx
import geonetworkx.settings as settings
import os
import shutil
from nose.tools import assert_is_instance, assert_equal, assert_in, assert_is
from nose.plugins.attrib import attr
import unittest
from shapely.geometry import Point
from itertools import chain


gnx_tu.SEED = 70595
NB_POINTS = 100


@attr('classes')
class TestClasses(unittest.TestCase):

    def setUp(self):
        file_dir = os.path.dirname(__file__)
        self.results_dir = os.path.join(file_dir, "datasets/results/")
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def tearDown(self):
        shutil.rmtree(self.results_dir)

    def test_default_node_geometry(self):
        """Try adding nodes without geometries and then change the geometry."""
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type):
                g = graph_type()
                g.add_nodes_from('Hello')
                for p in g.get_nodes_as_point_series():
                    assert_is(p, g.default_node_geometry)
                old_default_point = g.default_node_geometry
                new_default_point = Point(12, 21)
                g.default_node_geometry = new_default_point
                g.add_nodes_from([1, 2, 3])
                for n, p in g.get_nodes_as_point_series().items():
                    if str(n) in 'Hello':
                        assert_is(old_default_point, p)
                    else:
                        assert_is(new_default_point, p)

    def test_add_edge(self):
        """Test the add edge method for all graph types."""
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type):
                g = graph_type()
                g.add_edge(1, 2, geometry=gnx.LineString([(5, 4), (2, 7)]))
                assert_points_almost_equals(g.nodes[1]["geometry"], Point(5, 4))
                if g.is_multigraph():
                    new_key = g.add_edge(5, 3, 2, geometry=gnx.LineString([(5, 4), (9, 8)]))
                    assert_equal(new_key, 2, "Returned key is not the one provided")
                    assert_points_almost_equals(g.nodes[3]["geometry"], Point(9, 8))

    def test_add_edges_from(self):
        """Test the add edge method for all graph types."""
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type):
                g = graph_type()
                g.add_edges_from([(0, 1, dict(geometry=gnx.LineString([(0, 0), (1, 1)]))),
                                  (1, 2, dict(geometry=gnx.LineString([(1, 2), (2, 2)])))])
                assert_points_almost_equals(g.nodes[1]["geometry"], Point(1, 1))
                if g.is_multigraph():
                    g = graph_type()
                    new_keys = g.add_edges_from([(0, 1, 7, dict(geometry=gnx.LineString([(-1, 0), (1, 1)]))),
                                                 (1, 2, 8, dict(geometry=gnx.LineString([(1, 1), (2, 2)])))])
                    for k, vk in zip(new_keys, [7, 8]):
                        assert_equal(k, vk, "Returned key is not the one provided")
                    assert_points_almost_equals(g.nodes[2]["geometry"], Point(2, 2))

    def test_geograph_node_addition(self):
        """Try to add nodes an empty graph"""
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type):
                g = graph_type()
                g.add_nodes_from('Hello', **{g.nodes_geometry_key: Point(1, 2)})
                for n in 'Hello':
                    assert_in(n, g.nodes(), "Unable to find added node %s" % str(n))
                g.add_nodes_from([("H", {g.nodes_geometry_key: Point(2, 3)}),
                                  (2, {g.nodes_geometry_key: Point(5, 2)})])
                for n in chain('Hello', [2]):
                    assert_in(n, g.nodes(), "Unable to find added node %s" % str(n))
                assert_points_almost_equals(g.get_node_as_point("H"), Point(2, 3),
                                            "Replaced node geometry is not effective.")

    def test_graph_to_directed(self):
        graphs_directed_match = {gnx.GeoGraph: gnx.GeoDiGraph, gnx.GeoMultiGraph: gnx.GeoMultiDiGraph,
                                 gnx.GeoDiGraph: gnx.GeoDiGraph, gnx.GeoMultiDiGraph: gnx.GeoMultiDiGraph}
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                graph = get_random_geograph_subclass(NB_POINTS, graph_type)
                directed_graph = graph.to_directed()
                assert_is_instance(directed_graph, graphs_directed_match[graph_type])

    def test_graph_to_undirected(self):
        graphs_undirected_match = {gnx.GeoGraph: gnx.GeoGraph, gnx.GeoMultiGraph: gnx.GeoMultiGraph,
                                   gnx.GeoDiGraph: gnx.GeoGraph, gnx.GeoMultiDiGraph: gnx.GeoMultiGraph}
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                graph = get_random_geograph_subclass(NB_POINTS, graph_type)
                directed_graph = graph.to_undirected()
                assert_is_instance(directed_graph, graphs_undirected_match[graph_type])

    def test_crs_modification(self):
        graph = get_random_geograph_with_wgs84_scale(NB_POINTS, gnx.GeoMultiDiGraph)
        modified_graph = graph.to_crs(crs={'init': 'epsg:3857'}, inplace=False)
        re_modified_graph = modified_graph.to_crs(crs=gnx.settings.WGS84_CRS, inplace=False)
        assert_graphs_have_same_edges_geometry(graph, re_modified_graph, "Some edge geometries seems to be different"
                                                                         " after re-modification", tol=1e-2)
        assert_graphs_have_same_geonodes(graph, re_modified_graph, "Some nodes seems to be different after "
                                                                   "re-modification")

    def test_nodes_to_gdf(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_POINTS, graph_type)
                gdf = g.nodes_to_gdf()
                assert_equal(g.number_of_nodes(), len(gdf))

    def test_edges_to_gdf(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_POINTS, graph_type)
                gdf = g.edges_to_gdf()
                assert_equal(g.number_of_edges(), len(gdf))

    def test_spatial_keys_persistence(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_POINTS, graph_type)
                gnx.rename_edges_attribute(g, g.edges_geometry_key, "abcd")
                g.edges_geometry_key = "abcd"
                gnx.rename_nodes_attribute(g, g.nodes_geometry_key, "efgh")
                g.nodes_geometry_key = "efgh"
                g.crs = {'init': 'epsg:3945'}
                g2 = g.copy(as_view=False)
                assert_graphs_have_same_spatial_keys(g, g2)
                g3 = g.to_undirected(as_view=False)
                assert_graphs_have_same_spatial_keys(g, g3)
                g4 = g.to_directed(as_view=False)
                assert_graphs_have_same_spatial_keys(g, g4)

    def test_add_edges_gdf(self):
        for graph_type in ALL_CLASSES:
            with self.subTest(graph_type=graph_type, SEED=gnx_tu.SEED):
                g = get_random_geograph_subclass(NB_POINTS, graph_type)
                initial_edges = list(g.edges)
                nb_initial_edges = len(initial_edges)
                g2 = get_random_geograph_subclass(NB_POINTS, graph_type)
                edges_gdf = g2.edges_to_gdf()
                g.add_edges_from_gdf(edges_gdf, settings.EDGE_FIRST_NODE_COLUMN_NAME,
                                     gnx.settings.EDGE_SECOND_NODE_COLUMN_NAME)
                if g.is_multigraph():
                    g2_nb_edges = g2.number_of_edges()
                    assert_equal(g.number_of_edges(), g2_nb_edges + nb_initial_edges,
                                 "Resulting graph have not the right number of edges.")
                else:
                    g2_nb_distinct_edges = len([e for e in g2.edges if e not in initial_edges])
                    assert_equal(g.number_of_edges(), g2_nb_distinct_edges + nb_initial_edges,
                                 "Resulting graph have not the right number of edges.")
