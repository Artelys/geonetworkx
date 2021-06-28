# GeoNetworkX

Python tools for geographic graphs

[![CircleCI](https://circleci.com/gh/Artelys/geonetworkx.svg?style=svg)](https://circleci.com/gh/Artelys/geonetworkx)
[![Documentation Status](https://readthedocs.org/projects/geonetworkx/badge/?version=latest)](https://geonetworkx.readthedocs.io/en/latest/?badge=latest)
[![PyPI Package latest release](https://img.shields.io/pypi/v/geonetworkx.svg)](https://pypi.python.org/pypi/geonetworkx)

## Introduction

GeoNetworkX is a project to add support for geographic graphs to NetworkX (in the same way that GeoPandas support
geographic data to Pandas). It currently implements four data structures that extends the networkx graph classes (Graph,
MultiGraph, DiGraph, MultiDiGraph).


## Install

**Requirements**

* pyproj>=2.2
* geopy>=1.12.0
* geopandas>=0.7
* networkx>=2.3
* numpy>=1.12.0
* pandas>=0.25.0
* shapely>=1.2.18
* scipy>=0.19.0rc2
* nose>=1.3.7

Optional packages:

* srtm (for elevation data)
* pyvoronoi (for voronoi utils)
* osmnx (for OSM data)

**Installation**

```shell
    pip install geonetworkx
```

### Trouble when installing GeoNetworkX on Windows

If you are using `conda` on Windows, the binaries downloaded automatically
are broken, and `rtree` is
[unable to work properly](https://gis.stackexchange.com/questions/179706/installing-rtree-on-windows-64-bits).

A workaround is to download manually the binary from [this webpage](https://www.lfd.uci.edu/~gohlke/pythonlibs/#rtree).
Please download the binary corresponding to your system and your
Python version (3.6 or 3.7). You have notably to download:

- Binary for `Rtree`
- Binary for `GDAL`
- Binary for `Fiona`


## Documentation

Online documentation is available here: <https://geonetworkx.readthedocs.io>

## Tests

Tests can be launched with `unittest` with the following command:
```
python -m unittest discover -v geonetworkx
```
Or with `nose` like this:
```
nosetests geonetworkx -v --with-doctest
```

