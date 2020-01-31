
Supplement data
===============

Computing distances
-------------------

GeoNetworkX provides methods to compute distances within the coordinate
reference system of the graph. Typically, a method is given to add a length
attribute on edges. Different methods are available: euclidian distance but also
geodesic, great-circle or Vincenty distance (wrapped from `geopy`). All
available distances are stored within the dictionary `DISTANCE_MEASUREMENT_METHODS`.

.. doctest::

    >>> import geonetworkx as gnx
    >>> g = gnx.GeoGraph(crs=gnx.WGS84_CRS)
    >>> g.add_edge(1, 2, geometry=gnx.LineString([(-73.614, 45.504), (-73.632, 45.506)]))
    >>> gnx.fill_length_attribute(g)  # using geodesic distance
    >>> print(g.edges[(1, 2)]["length"])
    1424.174413518016
    >>> g.to_utm(inplace=True)
    >>> gnx.fill_length_attribute(g, only_missing=False)
    >>> print(g.edges[(1, 2)]["length"])  # using euclidian distance in UTM
    1423.8073619096585

A custom distance measurement method can be used by defining the appropriate
method in the settings. Here is an example implementing the Manhattan distance:

.. doctest::

    >>> def manhattan(p1, p2):
    ...     return abs(p1.x - p2.x) + abs(p1.y - p2.y)
    >>> gnx.settings.DISTANCE_MEASUREMENT_METHODS["manhattan"] = manhattan
    >>> gnx.fill_length_attribute(g, only_missing=False, method="manhattan")  # using manhattan distance
    >>> print(g.edges[(1, 2)]["length"])
    1608.0440213837428



Getting elevation data
----------------------

The elevation of nodes points can be filled as an attribute through the
`SRTM <https://en.wikipedia.org/wiki/Shuttle_Radar_Topography_Mission>`_ package.

.. doctest::

    >>> import geonetworkx as gnx
    >>> g = gnx.GeoGraph(crs=gnx.WGS84_CRS)
    >>> g.add_node(1, gnx.Point(5.145, 45.213))
    >>> gnx.fill_elevation_attribute(g)
    >>> print(g.nodes[1]["elevation[m]"])
    473


