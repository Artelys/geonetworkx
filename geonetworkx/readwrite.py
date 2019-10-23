# -*- coding: utf-8 -*-
import os
import pickle
import numpy as np
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point, LineString
from shapely.wkt import loads
from .geograph import GeoGraph
from .geomultigraph import GeoMultiGraph
from .geodigraph import GeoDiGraph
from .geomultidigraph import GeoMultiDiGraph
import geonetworkx.settings as settings
from geonetworkx.utils import get_crs_as_str


def parse_graph_as_geograph(graph, **attr):
    """Parse a ``networkx.Graph`` as a ``geonetworkx.GeoGraph`` with the closest geonetworkx graph type.

    Parameters
    ----------
    graph : nx.Graph, nx.DiGraph, nx.MultiGraph or nx.MultiDiGraph
        
    **attr :
        Potential spatial keys.

    Returns
    -------
    GeoGraph, GeoDiGraph, GeoMultiGraph or GeoMultiDiGraph
        Depending the orientation and multi edges properties.
    """
    if graph.is_directed():
        if graph.is_multigraph():
            geograph = GeoMultiDiGraph(graph, **attr)
        else:
            geograph = GeoDiGraph(graph, **attr)
    else:
        if graph.is_multigraph():
            geograph = GeoMultiGraph(graph, **attr)
        else:
            geograph = GeoGraph(graph, **attr)
    return geograph


def get_graph_with_wkt_geometry(geograph: GeoGraph) -> nx.Graph:
    """Modify the edges geometry attribute to a well-known text format to make the graph writable is some text formats.
    The returned graph is not as operational as the given one (edge geometries has been removed).

    Parameters
    ----------
    geograph: GeoGraph :
        Geograph to transform

    Returns
    -------
    nx.Graph
        A networkx graph with WKT geometries instead of shapely objects.

    See Also
    --------
    parse_nodes_attribute_as_wkt

    """
    graph_shallow_copy = geograph.to_nx_class()(geograph)
    nodes_geometries = nx.get_node_attributes(graph_shallow_copy, geograph.nodes_geometry_key)
    for n, point in nodes_geometries.items():
        if hasattr(point, "wkt"):
            graph_shallow_copy.nodes[n][geograph.nodes_geometry_key] = point.wkt
    edge_geometries = nx.get_edge_attributes(graph_shallow_copy, geograph.edges_geometry_key)
    for e, line in edge_geometries.items():
        if hasattr(line, "wkt"):
            graph_shallow_copy.edges[e][geograph.edges_geometry_key] = line.wkt
    return graph_shallow_copy


def parse_nodes_attribute_as_wkt(graph: nx.Graph, attribute_name: str):
    """Transform a graph nodes attribute from WKT to shapely objects. Attribute is replaced.

    Parameters
    ----------
    graph: nx.Graph :
        Graph to modify and parse
    attribute_name: str :
        Attribute to parse the nodes geometries

    See Also
    --------
    get_graph_with_wkt_geometry, parse_edges_attribute_as_wkt

    """
    wkt_points = nx.get_node_attributes(graph, attribute_name)
    for n, w in wkt_points.items():
        graph.nodes[n][attribute_name] = loads(w)


def parse_edges_attribute_as_wkt(graph: nx.Graph, attribute_name: str):
    """Transform a graph edges attribute from WKT to shapely objects. Attribute is replaced.

    Parameters
    ----------
    graph: nx.Graph :
        Graph to modify and parse
    attribute_name: str :
        Attribute to parse the edges geometries

    See Also
    --------
    get_graph_with_wkt_geometry, parse_nodes_attribute_as_wkt

    """
    wkt_lines = nx.get_edge_attributes(graph, attribute_name)
    for e, w in wkt_lines.items():
        graph.edges[e][attribute_name] = loads(w)


def stringify_crs(graph: GeoGraph):
    """Write the CRS attribute as a string."""
    if 'crs' in graph.graph and graph.graph['crs'] is not None:
        if not isinstance(graph.graph['crs'], str):
            graph.graph['crs'] = get_crs_as_str(graph.crs)


def read_gpickle(path, **attr):
    """Read geograph object in Python pickle format.

    Parameters
    ----------
    path : str
        Path where to read a graph object.
    **attr :
        Potential spatial keys.

    Returns
    -------
    GeoGraph, GeoDiGraph, GeoMultiGraph or GeoMultiDiGraph
        The parsed geograph.

    See Also
    --------
    write_gpickle, nx.read_gpickle, nx.write_gpickle

    """
    graph = nx.read_gpickle(path)
    if type(graph) in [nx.Graph, nx.DiGraph, nx.MultiDiGraph, nx.MultiGraph]:
        return parse_graph_as_geograph(graph, **attr)
    else:
        return graph


def write_gpickle(geograph, path, protocol=pickle.HIGHEST_PROTOCOL):
    """Write geograph object in Python pickle format.

    Parameters
    ----------
    geograph : GeoGraph, GeoDiGraph, GeoMultiGraph, GeoMultiDiGraph
        Geograph to write
    path :
        Path where to right the pickle file.
    protocol :
         See pickle protocols (Default value = pickle.HIGHEST_PROTOCOL).

    Returns
    -------

    See Also
    --------
    read_gpickle, nx.read_gpickle, nx.write_gpickle

    """
    nx.write_gpickle(geograph, path, protocol)


def read_graphml(path, node_type=str, edge_key_type=int, **attr) -> GeoGraph:
    """Read graph in GraphML format from path.

    Parameters
    ----------
    path :
        File path to the graphml file.
    node_type :
        See ``nx.read_graphml`` (Default value = str)
    edge_key_type :
        See ``nx.read_graphml`` (Default value = int)
    **attr :
        Potential spatial keys

    Returns
    -------
    GeoGraph, GeoDiGraph, GeoMultiGraph, GeoMultiDiGraph
        Parsed Geograph

    See Also
    --------
    write_graphml, nx.read_graphml, nx.write_graphml

    """
    graph = nx.read_graphml(path, node_type, edge_key_type)
    nodes_geometry_attr = attr.get("nodes_geometry_key",
                                   graph.graph.get("nodes_geometry_key",
                                                   settings.NODES_GEOMETRY_DEFAULT_KEY))
    parse_nodes_attribute_as_wkt(graph, nodes_geometry_attr)
    edges_geometry_attr = attr.get("edges_geometry_key",
                                   graph.graph.get("edges_geometry_key",
                                                   settings.EDGES_GEOMETRY_DEFAULT_KEY))
    parse_edges_attribute_as_wkt(graph, edges_geometry_attr)
    return parse_graph_as_geograph(graph, **attr)


def write_graphml(geograph, path, encoding='utf-8', prettyprint=True, infer_numeric_types=False):
    """Generate GraphML file for the given geograph.

    Parameters
    ----------
    geograph : GeoGraph, GeoDiGraph, GeoMultiGraph, GeoMultiDiGraph
        Geograph to write
    path : str
        Writing file path
    encoding :
         See ``nx.write_graphml`` (Default value = 'utf-8')
    prettyprint :
         See ``nx.write_graphml`` (Default value = True)
    infer_numeric_types :
         See ``nx.write_graphml`` (Default value = False)

    See Also
    --------
    read_graphml, nx.read_graphml, nx.write_graphml

    """
    graph_wkt = get_graph_with_wkt_geometry(geograph)
    stringify_crs(graph_wkt)
    nx.write_graphml(graph_wkt, path, encoding, prettyprint, infer_numeric_types)


def graph_nodes_to_gdf(graph: GeoGraph) -> gpd.GeoDataFrame:
    """Create and fill a GeoDataFrame (geopandas) from nodes of a networkX graph. The ``'geometry'`` attribute is used for
    shapes.

    Parameters
    ----------
    graph : GeoGraph
        Graph to parse

    Returns
    -------
    gpd.GeoDataFrame
        The resulting GeoDataFrame : one row is a node

    """
    nodes = {node: data for node, data in graph.nodes(data=True)}
    gdf_nodes = gpd.GeoDataFrame(nodes).T
    gdf_nodes['id'] = gdf_nodes.index
    gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
    if 'crs' in graph.graph:
        gdf_nodes.crs = graph.graph['crs']
    if 'name' in graph.graph:
        gdf_nodes.gdf_name = '{}_nodes'.format(graph.graph['name'])
    return gdf_nodes


def graph_edges_to_gdf(graph: nx.Graph) -> gpd.GeoDataFrame:
    """Create and fill a GeoDataFrame (geopandas) from edges of a networkX graph. The ``'geometry'`` attribute is used
     for shapes.

    Parameters
    ----------
    graph : nx.Graph
        Graph to parse

    Returns
    -------
    gpd.GeoDataFrame
        The resulting GeoDataFrame : one row is an edge

    """
    # create a list to hold our edges, then loop through each edge in the graph
    edges = []
    for u, v, data in graph.edges(data=True):
        # for each edge, add key and all attributes in data dict to the edge_details
        edge_details = {'u': u, 'v': v}
        for attr_key in data:
            edge_details[attr_key] = data[attr_key]
        # if edge doesn't already have a geometry attribute, create one now
        if 'geometry' not in data:
            point_u = Point((graph.nodes[u]['x'], graph.nodes[u]['y']))
            point_v = Point((graph.nodes[v]['x'], graph.nodes[v]['y']))
            edge_details['geometry'] = LineString([point_u, point_v])
        edges.append(edge_details)
    # create a GeoDataFrame from the list of edges and set the CRS
    gdf_edges = gpd.GeoDataFrame(edges)
    if 'crs' in graph.graph:
        gdf_edges.crs = graph.graph['crs']
    if 'name' in graph.graph:
        gdf_edges.gdf_name = '{}_edges'.format(graph.graph['name'])
    return gdf_edges


def parse_bool_columns_as_int(gdf: gpd.GeoDataFrame):
    """Transform bool columns into integer columns.

    Parameters
    ----------
    gdf: gpd.GeoDataFrame :
        GeoDataFrame to modify
    """
    for c in gdf.columns:
        if gdf[c].dtype == "bool":
            gdf[c] = gdf[c].astype("int")


def parse_numpy_types(gdf: gpd.GeoDataFrame):
    """Transform numpy types as scalar types.

    Parameters
    ----------
    gdf: gpd.GeoDataFrame :
        GeoDataFrame to modify
    """
    for c in gdf.columns:
        if any(map(lambda x: isinstance(x, np.generic), gdf[c])):
            for i in gdf.index:
                if isinstance(gdf.loc[i, c], np.generic):
                    gdf.loc[i, c] = np.asscalar(gdf.loc[i, c])


def stringify_unwritable_columns(gdf: gpd.GeoDataFrame):
    """Transform elements which have type bool or list to string

    Parameters
    ----------
    gdf: gpd.GeoDataFrame :
        GeoDataFrame to modify
    """
    valid_columns_types = ("int64", "float64")
    for c in gdf.columns:
        if not gdf[c].dtype in valid_columns_types and c != gdf._geometry_column_name:
            gdf[c] = list(map(str, gdf[c]))


def cast_for_fiona(gdf: gpd.GeoDataFrame):
    """Transform elements so that attributes can be writable by fiona.

    Parameters
    ----------
    gdf: gpd.GeoDataFrame :
        GeoDataFrame to modify
    """
    parse_bool_columns_as_int(gdf)
    parse_numpy_types(gdf)
    stringify_unwritable_columns(gdf)


def write_edges_to_geofile(graph: GeoGraph, file_name, driver="GPKG", fiona_cast=True):
    """Writes the edges of a geograph as a geographic file.

    Parameters
    ----------
    graph : GeoGraph, GeoDiGraph, GeoMultiGraph, GeoMultiDiGraph
        Graph to export
    file_name :
        File name (with path)
    driver :
        driver for export file format (GPKG, geojson, etc: can be found from ``fiona.supported_drivers``)
        (Default value = "GPKG")
    fiona_cast :
        If true, methods for casting types to writable fiona types are used (Default value = True)

    See Also
    --------
    write_geofile, write_nodes_to_geofile

    """
    gdf_edges = graph.edges_to_gdf()
    if fiona_cast:
        cast_for_fiona(gdf_edges)
    gdf_edges.to_file(file_name, driver=driver)


def write_nodes_to_geofile(graph: GeoGraph, file_name, driver="GPKG", fiona_cast=True):
    """Writes the nodes of a geograph as a geographic file.

    Parameters
    ----------
    graph : GeoGraph, GeoDiGraph, GeoMultiGraph, GeoMultiDiGraph
        Graph to export
    file_name :
        File name (with path)
    driver :
        driver for export file format (GPKG, geojson, etc: can be found from ``fiona.supported_drivers``)
        (Default value = "GPKG")
    fiona_cast :
        If true, methods for casting types to writable fiona types are used (Default value = True)

    See Also
    --------
    write_geofile, write_edges_to_geofile

    """
    gdf_nodes = graph.nodes_to_gdf()
    if fiona_cast:
        cast_for_fiona(gdf_nodes)
    gdf_nodes.to_file(file_name, driver=driver)


def write_geofile(graph: GeoGraph, path='./', nodes=True, edges=True, driver="GPKG", fiona_cast=True):
    """Export a networkx graph as a geographic files. Two files are generated: one for the nodes and one for the edges.
    The files names will be prefixed by the graph name and suffixed by "_edges" or "_nodes".

    Parameters
    ----------
    graph :
        Graph to export
    path :
        export directory (Default value = './')
    nodes :
        boolean to indicate whether export nodes or not. (Default value = True)
    edges :
        boolean to indicate whether export edges or not. (Default value = True)
    driver :
        driver for export file format (GPKG, geojson, etc: can be found from ``fiona.supported_drivers``)
         (Default value = "GPKG")
    fiona_cast :
        If true, methods for casting types to writable fiona types are used (Default value = True)

    See Also
    --------
    write_nodes_to_geofile, write_edges_to_geofile

    """
    if not os.path.exists(path):
        os.mkdir(path)
    if nodes:
        file_name = os.path.join(path, '{}_nodes'.format(graph.name))
        file_name += settings.KNOWN_FILES_EXTENSION.get(driver, "")
        write_nodes_to_geofile(graph, file_name, driver, fiona_cast)
    if edges:
        file_name = os.path.join(path, '{}_edges'.format(graph.name))
        file_name += settings.KNOWN_FILES_EXTENSION.get(driver, "")
        write_edges_to_geofile(graph, file_name, driver, fiona_cast)


def read_geograph_with_coordinates_attributes(graph: nx.Graph, x_key='x', y_key='y', **attr) -> GeoGraph:
    """Parse a `networkx` graph which have node's coordinates as attribute. This method can be useful to parse an output
    graph of the `osmnx` package.

    Parameters
    ----------
    graph : nx.Graph
        Given graph to parse. All nodes must have the ``x_key`` and ``y_key`` attributes.
    x_key :
        x-coordinates attribute to parse (Default value = 'x')
    y_key :
        y-coordinates attribute to parse (Default value = 'y')
    **attr :
        Optional geograph spatial keys.

    Returns
    -------
    GeoGraph, GeoDiGraph, GeoMultiGraph, GeoMultiDiGraph
        The parsed geograph (shallow copy of the input graph).

    """
    graph = graph.__class__(graph)
    x_coords = nx.get_node_attributes(graph, x_key)
    y_coords = nx.get_node_attributes(graph, y_key)
    nodes_geometry_key = attr.pop("nodes_geometry_key", settings.NODES_GEOMETRY_DEFAULT_KEY)
    edges_geometry_key = attr.pop("edges_geometry_key", settings.EDGES_GEOMETRY_DEFAULT_KEY)
    for n in graph.nodes:
        if n not in x_coords or n not in y_coords:
            raise ValueError("Unable to find coordinates for node : '%s'" % str(n))
        point = Point([x_coords[n], y_coords[n]])
        graph.nodes[n][nodes_geometry_key] = point
    graph.nodes_geometry_key = nodes_geometry_key
    graph.edges_geometry_key = edges_geometry_key
    return parse_graph_as_geograph(graph, **attr)


def read_geofiles(nodes_file_path: str, edges_file_path: str,
                  directed=True, multigraph=False,
                  node_index_attr=settings.NODE_ID_COLUMN_NAME,
                  edge_first_node_attr=settings.EDGE_FIRST_NODE_COLUMN_NAME,
                  edge_second_node_attr=settings.EDGE_SECOND_NODE_COLUMN_NAME):
    """Read geofiles to create a ``GeoGraph``. Geofiles are read with ``geopandas.read_file`` method.

    Parameters
    ----------
    nodes_file_path : str
        File path of nodes.
    edges_file_path : str
        File path of edges.
    directed : bool
        If ``True``, returns a directed graph. (Default value = True)
    multigraph : bool
        If ``True``, returns a multigraph. (Default value = False)
    node_index_attr : str
        Node id attribute in the geofile for nodes labeling. (Default value = settings.NODE_ID_COLUMN_NAME)
    edge_first_node_attr : str
        Edge first node attribute in the geofile. (Default value = settings.EDGE_FIRST_NODE_COLUMN_NAME)
    edge_second_node_attr : str
        Edge second node attribute in the geofile. (Default value = settings.EDGE_SECOND_NODE_COLUMN_NAME)

    Returns
    -------
    GeoGraph, GeoDiGraph, GeoMultiGraph, GeoMultiDiGraph
        A parsed ``Geograph``.

    See Also
    --------
    GeoGraph.add_nodes_from_gdf, GeoGraph.add_edges_from_gdf, geopandas.read_file

    """
    if directed:
        if multigraph:
            graph = GeoMultiDiGraph()
        else:
            graph = GeoDiGraph()
    else:
        if multigraph:
            graph = GeoMultiGraph()
        else:
            graph = GeoGraph()
    if nodes_file_path is not None:
        nodes_gdf = gpd.read_file(nodes_file_path)
        graph.nodes_geometry_key = nodes_gdf._geometry_column_name
        graph.crs = nodes_gdf.crs
        graph.add_nodes_from_gdf(nodes_gdf, node_index_attr)
    if edges_file_path is not None:
        edges_gdf = gpd.read_file(edges_file_path)
        if graph.crs is None and edges_gdf.crs is not None:
            graph.crs = edges_gdf.crs
        graph.add_edges_from_gdf(edges_gdf, edge_first_node_attr, edge_second_node_attr)
    return graph
