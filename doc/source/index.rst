.. GeoNetworkX documentation master file, created by
   sphinx-quickstart on Thu Jan 10 17:54:27 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


GeoNetworkX 0.5
===============

GeoNetworkX is a project to handle geospatial graphs. GeoNetworkX extends the NetworkX package to allow spatial
operations on geospatial graphs and benefit from the data structures and algorithm defined in NetworkX. Moreover, it
allows to use GeoPandas library tools on nodes and edges.

Description
-----------

The goal of GeoNetworkX is to embed a set of tools to handle geospatial graphs easily.
It combines the capabilities of networkx, geopandas and shapely, providing geospatial operations in networkx
high-level interface.


GeoNetworkX provides data structures that extends the networkx classes with this inheritance scheme:

.. figure:: ../figures/class_diagram.png
   :align: center

   GeoNetworkX inheritance graph

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   1_gettingStarted
   2_reading_and_writing_files
   3_supplement_data
   4_spatial_merge
   5_Isochrones
   reference/index


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



