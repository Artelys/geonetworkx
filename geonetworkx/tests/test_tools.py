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
os.chdir("geonetworkx/tests")
from nose.tools import assert_in


class TestTools():

    def test_spatial_points_merge(self):
        # test a spatial merge
        mdg = nx.read_gpickle("datasets/grenoble200_mdg.gpickle")
        graph = gnx.GeoMultiDiGraph(mdg)
        gnx.utils.fill_edges_missing_geometry_attributes(graph)
        points_gdf = gpd.read_file("datasets/grenoble200_buildings.geojson", driver="GeoJSON")
        gnx.tools.spatial_points_merge(graph, points_gdf)
        gnx.readwrite.export_graph_as_shape_file(graph, "datasets/results/")
        for p in points_gdf.index:
            assert_in(p, graph.nodes())
        nx.write_gpickle(graph, "datasets/grenoble100_merged.gpickle")
        gnx.utils.fill_length_attribute(graph)
        voro = nx.voronoi_cells(graph, set(points_gdf.index), weight='length')
        for ct, c in enumerate(voro):
            for n in voro[c]:
                graph.nodes[n]["voro"] = ct


    def test_spatial_graph_merge(self):
        streets_mdg = nx.read_gpickle("datasets/grenoble200_mdg.gpickle")
        streets_mdg = streets_mdg.to_undirected()
        base_graph = gnx.GeoMultiGraph(streets_mdg)
        gnx.utils.fill_edges_missing_geometry_attributes(base_graph)

        electrical_dg = nx.read_gpickle("datasets/grenoble200_electrical_dg.gpickle")
        electrical_dg = nx.MultiDiGraph(electrical_dg)
        original_edges = list(electrical_dg.edges(keys=True, data=True))
        electrical_mg = electrical_dg.to_undirected()
        other_graph = gnx.GeoMultiGraph(electrical_mg)
        gnx.order_well_lines(other_graph)
        gnx.join_lines_extremity_to_nodes_coordinates(other_graph)
        other_graph.graph["name"] = "electrical"
        gnx.export_graph_as_shape_file(other_graph, "datasets/results")
        merged_graph = gnx.tools.spatial_graph_merge(base_graph, other_graph, inplace=False)
        for n in other_graph.nodes():
            assert_in(n, merged_graph.nodes())
        for n in base_graph.nodes():
            assert_in(n, merged_graph.nodes())
        for e in other_graph.edges:
            assert_in(e, merged_graph.edges)


        merged_graph.graph["name"] = "merged_graph"
        gnx.readwrite.export_graph_as_shape_file(merged_graph, "datasets/results/")



