GeoNetworkX
===========

Python tools for geographic graphs


Introduction
------------

GeoNetworkX is a project to add support for geographic graphs to NetworkX (in the same way that GeoPandas support
geographic data to Pandas). It currently implements four data structures that extends the networkx graph classes (Graph,
 MultiGraph, DiGraph, MultiDiGraph).


Install
--------

**Requirements**

* pyproj>=1.9.6
* geopy>=1.12.0
* geopandas>=0.4.0
* networkx>=2.2
* numpy>=1.15.4
* shapely>=1.2.18
* scipy>=1.0.1

Optional packages:

* srtm (for elevation data)
* pyvoronoi (for voronoi utils)

**Installation**

``pip install geonetworkx``

Examples
--------

TODO

Tests
-----

Tests can be launched with `unittest` with the following command:
```
python -m unittest discover -v geonetworkx
```
Or with `nose` like this:
```
nosetests
```


