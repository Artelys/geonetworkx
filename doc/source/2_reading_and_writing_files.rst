
Reading and Writing Files
=========================

Reading Spatial Data
--------------------

As GeoNetworkX provides an interface to `geopandas` for the nodes and edges, it
is possible to read data from any vector-based spatial data supported by
`geopandas` (including ESRI Shapefile and GeoJSON).

Nodes and edges can be added to a given graph with the following methods:

.. code-block:: python

    import geonetworkx as gnx
    # Adding nodes and edges to an existing graph
    g = gnx.GeoGraph()
    g.add_nodes_from_gdf("copenhagen_streets_net_nodes.geojson")
    g.add_edges_from_gdf("copenhagen_streets_net_edges.geojson")
    gnx.read_geofiles(nodes_path, edges_path, directed=True, multigraph=True)

    # Creating a graph from existing files
    g = gnx.read_geofiles("copenhagen_streets_net_nodes.geojson",
                          "copenhagen_streets_net_edges.geojson",
                          edges_path, directed=True, multigraph=True)


Writing Spatial Data
--------------------

Geographs can be exported to same file formats as `geopandas`. Two files are
used to write a GeoGraph: one for nodes and one for edges. All the attributes of
the nodes and of the edges will be added in the files. If an attribute type is
not handled by `fiona` drivers, an attempt is made to cast it (see
`gnx.write_geofile` for more details).

.. code-block:: python

    g.name = "streets_graph"
    gnx.write_geofile(g, "test/path/", driver="GeoJSON")

The above code will write two GeoJSON files: :code:`test/path/streets_graph_nodes.geojson`
and :code:`test/path/streets_graph_edges.geojson` that can be directly read with
GIS software.


