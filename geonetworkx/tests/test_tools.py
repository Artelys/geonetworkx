# -*- coding: utf-8 -*-
import os
import shutil
import geopandas as gpd
import networkx as nx
import geonetworkx as gnx
import numpy as np
from shapely.geometry import Point
from nose.tools import assert_in, assert_less
import unittest
from nose.plugins.attrib import attr
import geonetworkx.testing.utils as gnx_tu
from geonetworkx.tests import datasets
import warnings

gnx_tu.SEED = 70595
data_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), "datasets")


@attr('tools')
class TestTools(unittest.TestCase):

    def setUp(self):
        file_dir = os.path.dirname(__file__)
        self.results_dir = os.path.join(file_dir, "datasets/results/")
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def tearDown(self):
        shutil.rmtree(self.results_dir)

    def test_spatial_points_merge(self):
        # test a spatial merge
        graph = datasets.get_copenhagen_street_net()
        gnx.utils.fill_edges_missing_geometry_attributes(graph)
        points_gdf = datasets.get_copenhagen_buildings()
        gnx.tools.spatial_points_merge(graph, points_gdf, inplace=True)
        for p in points_gdf.index:
            assert_in(p, graph.nodes())

    def test_spatial_graph_merge(self):
        base_graph = datasets.get_copenhagen_street_net().to_undirected()
        gnx.utils.fill_edges_missing_geometry_attributes(base_graph)
        base_graph.graph["name"] = "streets"

        other_graph = datasets.get_copenhagen_ferry_net().to_undirected()
        gnx.order_well_lines(other_graph)
        gnx.join_lines_extremity_to_nodes_coordinates(other_graph)
        other_graph.graph["name"] = "ferry"
        merged_graph = gnx.tools.spatial_graph_merge(base_graph, other_graph, inplace=False)
        merged_graph.graph["name"] = "merged_graph"
        for n in other_graph.nodes():
            assert_in(n, merged_graph.nodes())
        for n in base_graph.nodes():
            assert_in(n, merged_graph.nodes())
        for e in other_graph.edges:
            assert_in(e, merged_graph.edges)

    def test_spatial_graph_merge_with_non_distinct_nodes(self):
        # base graph definition
        nb_nodes = 20
        base_graph = gnx_tu.get_random_geograph(nb_nodes)
        gnx.fill_edges_missing_geometry_attributes(base_graph)

        # other graph definition
        nb_nodes_2 = 20
        edge_creation_prob = 0.1
        other_graph = gnx_tu.get_random_geograph(nb_nodes_2)
        nx.relabel_nodes(other_graph, {n: n + nb_nodes for n in other_graph.nodes}, False)
        distinct_nodes = list(other_graph.nodes())
        nx.set_node_attributes(other_graph, 2, "origin")
        common_nodes = list(range(nb_nodes // 4))
        common_part = nx.subgraph(base_graph, common_nodes)
        other_graph = nx.compose(other_graph, common_part)
        for n1 in common_nodes:
            for n2 in distinct_nodes:
                if np.random.rand() <= edge_creation_prob:
                    other_graph.add_edge(n1, n2)
        gnx.fill_edges_missing_geometry_attributes(other_graph)
        nx.set_node_attributes(other_graph, {n: 3 for n in common_nodes}, "origin")
        nx.set_node_attributes(base_graph, {n: 3 for n in common_nodes}, "origin")
        # nx.draw_networkx(base_graph, pos=nodes_coords, node_color ="red")
        # nx.draw_networkx(other_graph, pos=other_graph.get_nodes_coordinates(),
        #                  node_color=[2 if n in base_graph.nodes() else 3 for n in g2.nodes()])

        merged_graph = gnx.spatial_graph_merge(base_graph, other_graph, inplace=False)
        # origins = nx.get_node_attributes(merged_graph, "origin")
        # nc = [origins[n] if n in origins else 4 for n in merged_graph.nodes()]
        # nx.draw_networkx(merged_graph, pos=merged_graph.get_nodes_coordinates(), node_color=nc)

    def test_spatial_point_merge_with_filters(self):
        nb_nodes = 30
        graph = gnx_tu.get_random_geograph(nb_nodes)
        gnx.fill_edges_missing_geometry_attributes(graph)
        nb_points = 200
        points = np.random.rand(nb_points, 2)
        points_gdf = gpd.GeoDataFrame(columns=["geometry"])
        points_gdf["geometry"] = [Point([x, y]) for x, y in points]
        points_gdf.index = range(nb_nodes, nb_nodes + nb_points)
        # node filter
        node_filter = lambda n: n > 5
        merged_graph = gnx.spatial_points_merge(graph, points_gdf, node_filter=node_filter, inplace=False)
        for n in graph.nodes:
            assert_in(n, merged_graph.nodes)
        for n in points_gdf.index:
            assert_in(n, merged_graph.nodes)
        for n in graph.nodes:
            if not node_filter(n):
                edges = [e for e in graph.edges if n in e]
                for e in edges:
                    assert_in(e, merged_graph.edges)
        # edge filter
        edge_filter = lambda u, v: u not in range(5, 10) and v not in range(12, 20)
        merged_graph = gnx.spatial_points_merge(graph, points_gdf, edge_filter=edge_filter, inplace=False)
        for n in graph.nodes:
            assert_in(n, merged_graph.nodes)
        for n in points_gdf.index:
            assert_in(n, merged_graph.nodes)
        for e in graph.edges:
            if not edge_filter(e[0], e[1]):
                assert_in(e, merged_graph.edges)
        # node filter and edge filter
        merged_graph = gnx.spatial_points_merge(graph, points_gdf, edge_filter=edge_filter, node_filter=node_filter,
                                                inplace=False)
        for n in graph.nodes:
            assert_in(n, merged_graph.nodes)
        for n in points_gdf.index:
            assert_in(n, merged_graph.nodes)
        for e in graph.edges:
            if not edge_filter(e[0], e[1]):
                assert_in(e, merged_graph.edges)
        for n in graph.nodes:
            if not node_filter(n):
                edges = [e for e in graph.edges if n in e]
                for e in edges:
                    assert_in(e, merged_graph.edges)

    def test_isochrone(self):
        if not gnx_tu.check_optional_package_presence("pyvoronoi"):
            warnings.warn("Missing optional package for test 'test_isochrone': pyvoronoi")
            return
        # Read data
        gmg = datasets.get_grenoble_streets_500()
        gnx.fill_edges_missing_geometry_attributes(gmg)
        gnx.fill_length_attribute(gmg, "length", True)

        # Compute Isochrone polygon
        source = 312173744
        limit = 400  # meters
        iso_polygon = gnx.isochrone_polygon(gmg, source, limit, "length")
        # Try to compute the distance for a set of point that are in the polygon bounding box
        nb_test_points = 20
        np.random.seed(gnx_tu.SEED)
        bb = iso_polygon.bounds
        points_coords = np.random.rand(nb_test_points, 2) * np.array([bb[2] - bb[0], bb[3] - bb[1]]) + \
                        np.array([bb[0], bb[1]])
        points = gpd.GeoDataFrame([Point(p) for p in points_coords], columns=["geometry"], crs=gmg.crs)
        merged_graph = gnx.spatial_points_merge(gmg, points_gdf=points, inplace=False)
        gnx.fill_length_attribute(merged_graph, "length", True)
        tol = 0.05 * limit
        for p in points.index:
            shortest_distance = nx.shortest_path_length(merged_graph, source, p, weight="length")
            if points.at[p, "geometry"].intersects(iso_polygon):
                assert_less(shortest_distance - tol, limit, "The shortest path between a point intersecting the iso-polygon"
                                                            " must be shorter than the isochrone limit.")
            else:
                assert_less(limit, shortest_distance + tol, "The shortest path between a point not intersecting the "
                                                            "iso-polygon must be longer than the isochrone limit.")
