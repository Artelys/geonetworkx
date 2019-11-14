
Getting started
===============

Installation
------------

GeoNetworkX can be installed with pip with the following command:

.. code-block:: shell

    pip install geonetworkx

.. warning::
    GeoNetworkX needs packages that have C dependencies that may need to compiled and installed manually
    (`shapely`, `fiona`, `pyproj` and `rtree`). For Windows users, wheels can be found at `Christopher Gohlke's
    website <https://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.

What's a GeoGraph ?
-------------------

A geograph is an object extending classical graph definition with topological space. For example, it can be a road
network where edges represent streets and nodes represent their intersection. Other applications can be found in
electrical networks, railway networks, etc.
Mathematically, a geograph is defined with the following elements:

    * A topological space :math:`S` with a distance measurement application :math:`d: S \times S \rightarrow \mathbb{R}^+`
    * A graph :math:`G(N, E)` with :math:`N` a finite set of vertices and :math:`E \subset N^2` at set of pairs of vertices.
    * :math:`P := \bigcup_{n \in N} p_n` with :math:`p_n \in S` the coordinates of the node :math:`n`.
    * :math:`L := \bigcup_{(u, v) \in E} l_{u, v}` with :math:`l_{u, v} \subset S` a topological curve starting at :math:`p_u` and ending at :math:`p_v`.

The space :math:`S` is usually here considered to be :math:`\mathbb{R}^2` with the euclidian distance, or the WGS84
spheroid with the great-circle distance (or Vincenty distance).

Closest edge rule
-----------------

The implementation uses the closest edge rule to connect a topological point :math:`p \in S` to a geograph. This rule
define a connection point :math:`i_p` :

    .. math::
        i_p := \text{proj}_{L}(p) = \text{argmin}\{d(p, x) | x \in L\}

This rule allows to connect any point of the topological space to the geograph. In the street network example, it means
finding the closest street for starting a trip.

Implementation details
----------------------

The implementation of geographs in GeoNetworkX is based on the following hypothesis:

    * All nodes have coordinates stored in a ``shapely.geometry.Point`` object.
    * Edges may have geometry stored in a ``shapely.geometry.LineString`` object.
    * A geograph may have a Coordinate Reference System (CRS) using GeoPandas implementation.

An edge may not have a geometry but it is supposed that it can be deduced by a simple "straight" line between the two
nodes.