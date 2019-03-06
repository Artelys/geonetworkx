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
    """Parse a ``networkx.Graph`` as a ``geonetworkx.GeoGraph`` with the closest geonetworkx graph type."""
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
    The returned graph is not as operational as the given one (edge geometries has been removed)"""
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
    """Parse nodes geometries with a wkt (well known text) attribute"""
    wkt_points = nx.get_node_attributes(graph, attribute_name)
    for n, w in wkt_points.items():
        graph.nodes[n][attribute_name] = loads(w)


def parse_edges_attribute_as_wkt(graph: nx.Graph, attribute_name: str):
    """Parse edge geometries with a wkt (well known text) attribute."""
    wkt_lines = nx.get_edge_attributes(graph, attribute_name)
    for e, w in wkt_lines.items():
        graph.edges[e][attribute_name] = loads(w)


def stringify_crs(graph: GeoGraph):
    """Write the CRS attribute as a string."""
    if 'crs' in graph.graph and graph.graph['crs'] is not None:
        if not isinstance(graph.graph['crs'], str):
            graph.graph['crs'] = get_crs_as_str(graph.crs)


def read_gpickle(path, **attr):
    """Read graph object in Python pickle format."""
    graph = nx.read_gpickle(path)
    if isinstance(graph, (nx.Graph, nx.DiGraph, nx.MultiDiGraph, nx.MultiGraph)):
        return parse_graph_as_geograph(graph, **attr)
    else:
        return graph


def write_gpickle(geograph, path, protocol=pickle.HIGHEST_PROTOCOL):
    """Write graph object in Python pickle format."""
    nx.write_gpickle(geograph, path, protocol)


def read_graphml(path, node_type=str, edge_key_type=int, **attr):
    """Read graph in GraphML format from path.

    :param path: File path to the graphml file.
    :param node_type: See ``nx.read_graphml``
    :param edge_key_type: See ``nx.read_graphml``
    :param attr: GeoGraph spatial keys specification
    :return: 
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
    """Generate GraphML lines for G"""
    graph_wkt = get_graph_with_wkt_geometry(geograph)
    stringify_crs(graph_wkt)
    nx.write_graphml(graph_wkt, path, encoding, prettyprint, infer_numeric_types)


def graph_nodes_to_gdf(graph: GeoGraph) -> gpd.GeoDataFrame:
    """
    Create and fill a GeoDataFrame (geopandas) from nodes of a networkX graph. The ``'geometry'`` attribute is used for
    shapes.

    :param graph: Graph to parse
    :return: The resulting GeoDataFrame : one row is a node
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
    """
    Create and fill a GeoDataFrame (geopandas) from edges of a networkX graph. The ``'geometry'`` attribute is used for
    shapes.

    :param graph: Graph to parse
    :return: The resulting GeoDataFrame : one row is an edge
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
    """Transform bool columns into integer columns."""
    for c in gdf.columns:
        if gdf[c].dtype == "bool":
            gdf[c] = gdf[c].astype("int")


def parse_numpy_types(gdf: gpd.GeoDataFrame):
    """Transform numpy types as scalar types."""
    for c in gdf.columns:
        if any(map(lambda x: isinstance(x, np.generic), gdf[c])):
            for i in gdf.index:
                if isinstance(gdf.loc[i, c], np.generic):
                    gdf.loc[i, c] = np.asscalar(gdf.loc[i, c])


def stringify_unwritable_columns(gdf: gpd.GeoDataFrame):
    """Transform elements which have type bool or list to string"""
    valid_columns_types = ("int64", "float64")
    for c in gdf.columns:
        if not gdf[c].dtype in valid_columns_types and c != settings.GPD_GEOMETRY_KEY:
            gdf[c] = list(map(str, gdf[c]))


def cast_for_fiona(gdf: gpd.GeoDataFrame):
    """Transform elements so that attributes can be writable by fiona."""
    parse_bool_columns_as_int(gdf)
    parse_numpy_types(gdf)
    stringify_unwritable_columns(gdf)


def write_edges_to_geofile(graph: GeoGraph, file_name, driver="GPKG", fiona_cast=True):
    """ Writes the edges of a geograph as a geographic file.

    :param graph: Graph to export
    :param file_name: File name (with path)
    :param driver: driver for export file format (GPKG, geojson, etc: can be found from ``fiona.supported_drivers``)
    :param fiona_cast: If true, methods for casting types to writable fiona types are used
    """
    gdf_edges = graph.edges_to_gdf()
    if fiona_cast:
        cast_for_fiona(gdf_edges)
    gdf_edges.to_file(file_name, driver=driver)


def write_nodes_to_geofile(graph: GeoGraph, file_name, driver="GPKG", fiona_cast=True):
    """ Writes the nodes of a geograph as a geographic file.

    :param graph: Graph to export
    :param file_name: File name (with path)
    :param driver: driver for export file format (GPKG, geojson, etc: can be found from ``fiona.supported_drivers``)
    :param fiona_cast: If true, methods for casting types to writable fiona types are used
    """
    gdf_nodes = graph.nodes_to_gdf()
    if fiona_cast:
        cast_for_fiona(gdf_nodes)
    gdf_nodes.to_file(file_name, driver=driver)


def write_geofile(graph: GeoGraph, path='./', nodes=True, edges=True, driver="GPKG", fiona_cast=False):
    """
    Export a networkx graph as a geographic file.

    :param graph: Graph to export
    :param path: export directory
    :param nodes: boolean to indicate whether export nodes or not.
    :param edges: boolean to indicate whether export edges or not.
    :param driver: driver for export file format (GPKG, geojson, etc: can be found from ``fiona.supported_drivers``)
    :param fiona_cast: If true, methods for casting types to writable fiona types are used
    :return: None
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
    """Parse a `networkx` graph which have node's coordinates as attribute.

    :param graph: Given graph to parse. All nodes must have the ``x_key`` and ``y_key`` attributes.
    :param x_key: x-coordinates attribute to parse
    :param y_key: y-coordinates attribute to parse
    :param attr: Optional geograph spatial keys.
    :return: The parsed geograph (shallow copy of the input graph).
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

