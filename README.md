# GeoNetworkX

Python tools for geographic graphs


## Introduction

GeoNetworkX is a project to add support for geographic graphs to NetworkX (in the same way that GeoPandas support
geographic data to Pandas). It currently implements four data structures that extends the networkx graph classes (Graph,
MultiGraph, DiGraph, MultiDiGraph).


## Install

**Requirements**

* Fiona==1.8.6
* geopy==1.2.0
* geopandas==0.4.0
* networkx==2.3
* numpy>=1.15.4
* pyproj==1.9.6
* rtree==0.8.3
* scipy>=1.0.1
* shapely>=1.2.18

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
- Binary for `Fiona`
- Binary for `PyProj 1.9.6` (mandatory version!)



## Examples

TODO

## Tests

Tests can be launched with `unittest` with the following command:
```
python -m unittest discover -v geonetworkx
```
Or with `nose` like this:
```
nosetests
```

