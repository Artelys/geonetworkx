import networkx as nx

def remove_isolates(graph: nx.Graph) -> int:
    """
    Removes all isolates nodes in the given graph.
    :param graph: A graph on which to remove all isolates
    :return: Number of removed isolates
    """
    isolates = list(nx.isolates(graph))
    graph.remove_nodes_from(isolates)
    return len(isolates)


def remove_self_loop_edges(graph: nx.Graph) -> int:
    """
    Remove self loop edges on nodes of the given graph.
    :param graph: A graph on which to remove all self loops.
    :return: The number of removed self loops
    """
    self_loops_edges = list(nx.selfloop_edges(graph))
    graph.remove_edges_from(self_loops_edges)
    return len(self_loops_edges)


def remove_small_connected_components(graph: nx.Graph, minimum_allowed_size: int) -> int:
    """
    Remove all connected components having strictly less than 'minimum_allowed_size'
    :param graph: The graph on which to remove connected components
    :param minimum_allowed_size: The minimum number of nodes where a connected component is kept.
    :return: The number of removed connected components
    """
    connected_components = list(nx.connected_components(graph))
    nb_removed_cc = 0
    for c_ix, cc in enumerate(connected_components):
        if len(cc) < minimum_allowed_size:
            graph.remove_nodes_from(cc)
            nb_removed_cc += 1
        else:
            for n in cc:
                graph.nodes[n]['cc'] = c_ix
    return nb_removed_cc
