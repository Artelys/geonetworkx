
Isochrones
==========

Isochrones defines the reach of a location within a distance limit. An isochrone polygon is defined by all reachable
points from a source node trough a given geograph.

Using an additional distance function :math:`d_G: N \times N \rightarrow \mathbb{R}^+` corresponding to a shortest path
distance between two nodes in the graph, the isochrone polygon can be defined as:

    .. math::
        I_n^\epsilon := \{x \in S: d(i_x, x) + d_G(n, i_x) \leq \epsilon\}

with :math:`i_x := \text{proj}_L(x)`

# TODO: simplify math


The core method is based on Shortest Path Tree generation (SPT) (or ego-graph as in NetworkX). This tree contains all
nodes reachable within the :math:`\epsilon` distance from the source node. To get an isochrone polygon approximation,
this tree has to be "buffered" to represent the boundaries of the SPT.
Such polygons can be approximate by various methods (see `Isochrones OSM wiki <https://wiki.openstreetmap.org/wiki/Isochrone>`_).
Two methods have been implemented in GeoNetworkX:

    * :math:`\alpha`-shapes for fast approximation
    * Natural neighborhood with Voronoi edges cells computation for a sharp approximation


Voronoi edges cells
-------------------

This methods tries to find


:math:`\alpha`-shape
--------------------


Test