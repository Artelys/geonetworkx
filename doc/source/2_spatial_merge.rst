
Spatial merge tools
===================


GeoNetworkX implements methods for map-matching points to the geograph. That is to say finding, for a query point, the
closest edge or node in a geograph. Mathematically, it means solving the following optimization problem for a query
point :math:`p \in S`:

    .. math::
        \min_{x \in L} d(p, x)

A frequently encountered problem is finding the closest edge in a geograph for a set of points :math:`P`. If done
naively, a nested loop on points and edges is performed to compute all distances. This introduces a high computational
cost (:math:`o(|P|\times|E|)`) that can be avoided by using the right data structure. GeoNetworkX uses kd-trees to
efficiently solve this problem. Theses allow to find the optimal solution without having to compute all distances
(:math:`o((|P| + |E|) \log|E| )`).


To compute the closest edge, all edge geometries are discretized within a tolerance distance :math:`\epsilon > 0`. This
method is not exact if the coordinate of the geograph are not unprojected (using latitude and longitude angles), but
produces fairly good results. The implemented method uses kd-trees implemented in the
`Scipy Spatial <https://docs.scipy.org/doc/scipy/reference/spatial.html>`_ package (see
`cKDTree <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html#scipy.spatial.cKDTree>`_).


Spatial points merge
--------------------

GeoNetworkX implements a method to add a set of point to a geograph as nodes using the closest edge rule
(``gnx.spatial_points_merge``). This method not only find the closest edge but generate the new edges to connect the
new nodes to the geograph:

.. figure:: ../figures/building_projection_graph2.png
    :align: center
    :figclass: align-center

    Illustration of work done in the ``spatial_points_merge`` method. Initial geograph edges are in black, nodes in blue
    , new nodes are in green, new intersection nodes in red.


Spatial graph merge
-------------------

An additional useful feature that provides GeoNetworkX is geographs merge. That is to say, from a base graph, adding
another graph on top of it and setting the right edges connect both.
This feature may be very useful for multimodal transport routing. For example, a use case is to merge a street graph
with a subway system graph to find an optimal route combining walk and subway transportation. To do so, the closest
street of each subway station has to be found and an edge has to be added to link them. This is what is done in the
``gnx.spatial_graph_merge`` method.

Here is a use case using this method:

```

```


