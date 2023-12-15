import geopandas as gpd
import networkx as nx
from shapely import set_precision
from shapely.geometry import Point


def lines_to_graph(lines_gdf):
    """
    Converts a LineString gpd.GeoDataFrame to a nx.MultiDiGraph.
    Columns are preserved as edge attributes, node labels are
    initially set to the respective coordinate tuples.
    inspired by momepy package > utils.gdf_to_nx() function

    Parameters:
            lines_gdf (gpd.GeoDataFrame): LineStrings plus additional metadata

    Returns:
            graph (nx.MultiDiGraph): graph representation of lines_gdf
    """

    network_gdf = lines_gdf.copy()

    graph = nx.MultiDiGraph()
    graph.graph["crs"] = network_gdf.crs

    network_gdf["length"] = network_gdf.geometry.length
    if len(network_gdf.geom_type.unique()) > 1:
        network_gdf = network_gdf.explode(index_parts=True)

    network_gdf["geometry"] = set_precision(network_gdf.geometry, 0.1)
    cols = list(network_gdf.columns)

    key = 0
    vertices = []
    for row in network_gdf.itertuples():
        first = row.geometry.coords[0]
        last = row.geometry.coords[-1]

        data = list(row)[1:]
        attributes = dict(zip(cols, data))

        if (first, last) in vertices:
            key += 1
        graph.add_edge(first, last, key=key, **attributes)
        vertices.append((first, last))

    nodes = set(sum(vertices, ()))
    for node in nodes:
        graph.nodes[node]["geometry"] = Point(node)

    return graph


def graph_to_gdfs(graph):
    """
    Convert a nx.MultiDiGraph to node and/or edge gpd.GeoDataFrame.
    inspired by osmnx package > utils_graph.graph_to_gdfs() function

    Parameters:
            graph (nx.MultiDiGraph): graph representation of lines_gdf

    Returns:
            gpd.GeoDataFrames or tuple: gdf_nodes or gdf_edges or tuple
            of (gdf_nodes, gdf_edges). gdf_nodes is indexed by (x, y)
            coordinate and gdf_edges is multi-indexed by (u, v, key)
            following normal nx.MultiDiGraph structure.
    """

    crs = graph.graph["crs"]

    if not graph.nodes:
        raise ValueError("Graph contains no nodes.")

    nodes, data = zip(*graph.nodes(data=True))
    geom = list(d["geometry"] for d in data)
    gdf_nodes = gpd.GeoDataFrame(data, index=nodes, crs=crs, geometry=geom)

    if not graph.edges:
        raise ValueError("Graph contains no edges.")

    u, v, k, data = zip(*graph.edges(keys=True, data=True))
    geom = list(d["geometry"] for d in data)
    gdf_edges = gpd.GeoDataFrame(data, crs=crs, geometry=geom)

    gdf_edges["u"] = u
    gdf_edges["v"] = v
    gdf_edges["k"] = k
    gdf_edges.set_index(["u", "v", "k"], inplace=True)

    return gdf_nodes, gdf_edges
