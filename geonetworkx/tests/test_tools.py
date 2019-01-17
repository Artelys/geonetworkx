"""
    File name: test_tools
    Author: Artelys
    Creation date: 09/01/2019
    Python Version: 3.6
"""
import sys, os
import geopandas as gpd
import networkx as nx
import geonetworkx as gnx
import numpy as np
from shapely.geometry import Point
#os.chdir("geonetworkx/tests")
from nose.tools import assert_in, assert_true

SEED = 70595
np.random.seed(SEED)
data_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), "datasets")


class TestTools():
    def test_spatial_points_merge(self):
        # test a spatial merge
        mdg = nx.read_gpickle(os.path.join(data_directory, "grenoble200_mdg.gpickle"))
        graph = gnx.GeoMultiDiGraph(mdg)
        gnx.utils.fill_edges_missing_geometry_attributes(graph)
        points_gdf = gpd.read_file(os.path.join(data_directory, "grenoble200_buildings.geojson"), driver="GeoJSON")
        gnx.tools.spatial_points_merge(graph, points_gdf)
        #gnx.readwrite.export_graph_as_shape_file(graph, "datasets/results/")
        for p in points_gdf.index:
            assert_in(p, graph.nodes())
        #nx.write_gpickle(graph, "datasets/grenoble100_merged.gpickle")
        gnx.utils.fill_length_attribute(graph)
        voro = nx.voronoi_cells(graph, set(points_gdf.index), weight='length')
        for ct, c in enumerate(voro):
            for n in voro[c]:
                graph.nodes[n]["voro"] = ct


    def test_spatial_graph_merge(self):
        streets_mdg = nx.read_gpickle(os.path.join(data_directory, "grenoble200_mdg.gpickle"))
        streets_mdg = streets_mdg.to_undirected()
        base_graph = gnx.GeoMultiGraph(streets_mdg)
        gnx.utils.fill_edges_missing_geometry_attributes(base_graph)
        base_graph.graph["name"] = "streets"
        #gnx.export_graph_as_shape_file(base_graph, "datasets/results/")

        electrical_dg = nx.read_gpickle(os.path.join(data_directory, "grenoble200_electrical_dg.gpickle"))
        electrical_dg = nx.MultiDiGraph(electrical_dg)
        original_edges = list(electrical_dg.edges(keys=True, data=True))
        electrical_mg = electrical_dg.to_undirected()
        other_graph = gnx.GeoMultiGraph(electrical_mg)
        gnx.order_well_lines(other_graph)
        gnx.join_lines_extremity_to_nodes_coordinates(other_graph)
        other_graph.graph["name"] = "electrical"
        #gnx.export_graph_as_shape_file(other_graph, "datasets/results")
        merged_graph = gnx.tools.spatial_graph_merge(base_graph, other_graph, inplace=False)
        merged_graph.graph["name"] = "merged_graph"
        #gnx.readwrite.export_graph_as_shape_file(merged_graph, "datasets/results/")
        for n in other_graph.nodes():
            assert_in(n, merged_graph.nodes())
        for n in base_graph.nodes():
            assert_in(n, merged_graph.nodes())
        for e in other_graph.edges:
            assert_in(e, merged_graph.edges)



    def test_spatial_graph_merge_with_non_distinct_nodes(self):
        # base graph definition
        nb_nodes = 20
        edge_creation_prob = 0.1
        g = nx.fast_gnp_random_graph(nb_nodes, edge_creation_prob, seed=SEED, directed=False)
        nodes_coords = nx.kamada_kawai_layout(g)
        nx.set_node_attributes(g, {n: coords[0] for n, coords in nodes_coords.items()}, 'x')
        nx.set_node_attributes(g, {n: coords[1] for n, coords in nodes_coords.items()}, 'y')
        nx.set_node_attributes(g, 1, "origin")
        base_graph = gnx.GeoGraph(g)
        gnx.fill_edges_missing_geometry_attributes(base_graph)

        # other graph definition
        nb_nodes_2 = 20
        g2 = nx.fast_gnp_random_graph(nb_nodes_2,edge_creation_prob, seed=SEED + 1, directed=False)
        nx.relabel_nodes(g2, {n: n + nb_nodes for n in g2.nodes}, False)
        distinct_nodes = list(g2.nodes())
        nodes_coords_2 = nx.kamada_kawai_layout(g2)
        nx.set_node_attributes(g2, {n: coords[0] for n, coords in nodes_coords_2.items()}, 'x')
        nx.set_node_attributes(g2, {n: coords[1] for n, coords in nodes_coords_2.items()}, 'y')
        nx.set_node_attributes(g2, 2, "origin")
        common_nodes = list(range(nb_nodes // 4))
        common_part = nx.subgraph(base_graph, common_nodes)
        g2 = nx.compose(g2, common_part)
        for n1 in common_nodes:
            for n2 in distinct_nodes:
                if np.random.rand() <= edge_creation_prob:
                    g2.add_edge(n1, n2)
        other_graph = gnx.GeoGraph(g2)
        gnx.fill_edges_missing_geometry_attributes(other_graph)
        nx.set_node_attributes(other_graph, {n: 3 for n in common_nodes}, "origin")
        nx.set_node_attributes(base_graph, {n: 3 for n in common_nodes}, "origin")
        nx.draw_networkx(base_graph, pos=nodes_coords, node_color ="red")
        nx.draw_networkx(other_graph, pos=other_graph.get_nodes_coordinates(), node_color=[2 if n in base_graph.nodes() else 3 for n in g2.nodes()])

        merged_graph = gnx.spatial_graph_merge(base_graph, other_graph, inplace=False)
        origins = nx.get_node_attributes(merged_graph, "origin")
        nc = [origins[n] if n in origins else 4 for n in merged_graph.nodes()]
        nx.draw_networkx(merged_graph, pos=merged_graph.get_nodes_coordinates(), node_color=nc)

    def test_spatial_point_merge_with_filters(self):
        nb_nodes = 30
        edge_creation_prob = 0.1
        g = nx.fast_gnp_random_graph(nb_nodes, edge_creation_prob, seed=SEED, directed=False)
        nodes_coords = nx.kamada_kawai_layout(g)
        nx.set_node_attributes(g, {n: coords[0] for n, coords in nodes_coords.items()}, 'x')
        nx.set_node_attributes(g, {n: coords[1] for n, coords in nodes_coords.items()}, 'y')
        graph = gnx.GeoGraph(g)
        gnx.fill_edges_missing_geometry_attributes(graph)
        nb_points = 200
        points = np.random.rand(nb_points, 2)
        points_gdf = gpd.GeoDataFrame(columns=["geometry"])
        points_gdf["geometry"] = [Point([x,y]) for x,y in points]
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
        merged_graph = gnx.spatial_points_merge(graph, points_gdf, edge_filter=edge_filter, node_filter=node_filter, inplace=False)
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


