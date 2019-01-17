"""
    File name: test_classes
    Author: Artelys
    Creation date: 17/01/2019
    Python Version: 3.6
"""
from geonetworkx.testing import get_random_geograph, get_random_geomultigraph
import geonetworkx as gnx
import numpy as np
from nose.tools import assert_is_instance


SEED = 70595
np.random.seed(SEED)


class TestClasses():
    def test_graph_to_directed(self):
        nb_nodes = 200
        graph = get_random_geograph(nb_nodes, SEED)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoDiGraph)

    def test_mutligraph_to_directed(self):
        nb_nodes = 200
        graph = get_random_geomultigraph(nb_nodes, SEED + 1)
        directed_graph = graph.to_directed()
        assert_is_instance(directed_graph, gnx.GeoMultiDiGraph)

