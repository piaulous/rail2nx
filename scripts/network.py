import networkx as nx
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
