import logging
from itertools import product

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely import line_merge, shortest_line
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import nearest_points, snap, split

from scripts.tools import geom_to_int, read_config, tuple_coords

cfg = read_config()
logger = logging.getLogger("rail2nx")
pd.options.mode.chained_assignment = None


def lines_to_graph(lines_gdf, line_repair=True, keep_isolates=False):
    """
    Converts a LineString gpd.GeoDataFrame to a nx.Graph.
    Columns are preserved as edge attributes, node labels are
    initially set to the respective coordinate tuples.
    inspired by momepy package > utils.gdf_to_nx() function

    Parameters:
            lines_gdf (gpd.GeoDataFrame): LineStrings plus additional metadata

    Returns:
            graph (nx.Graph): graph representation of lines_gdf
    """

    network_gdf = lines_gdf.copy()

    if line_repair:
        network_gdf = repair_lines(network_gdf)

    graph = nx.Graph()
    graph.graph["crs"] = network_gdf.crs

    network_gdf["length"] = network_gdf.geometry.length
    if len(network_gdf.geom_type.unique()) > 1:
        network_gdf = network_gdf.explode(index_parts=True)
        network_gdf["length"] = network_gdf.geometry.length

    # network_gdf = round_gdf(network_gdf)
    cols = list(network_gdf.columns)

    vertices = []
    for row in network_gdf.itertuples():
        first = tuple_coords(row.geometry.coords[0])
        last = tuple_coords(row.geometry.coords[-1])

        data = list(row)[1:]
        attributes = dict(zip(cols, data))

        graph.add_edge(first, last, **attributes)
        vertices.append((first, last))

    nodes = set(sum(vertices, ()))
    for node in nodes:
        graph.nodes[node]["geometry"] = Point(node)

    if not keep_isolates:
        graph = remove_isolates(graph)

    return graph


def repair_lines(gdf_lines):
    line_gaps = pd.read_csv("data/line_gaps.csv")
    for u, v in zip(line_gaps.FID_rail_d_1, line_gaps.FID_rail_d_2):
        geoms_1 = gdf_lines[gdf_lines.FID_rail_d == u]
        geoms_2 = gdf_lines[gdf_lines.FID_rail_d == v]
        if any([geoms_1.empty, geoms_2.empty]):
            continue

        line_dict = {}
        for idx_1, idx_2 in product(geoms_1.index, geoms_2.index):
            line_1 = gdf_lines.loc[idx_1, "geometry"]
            line_2 = gdf_lines.loc[idx_2, "geometry"]
            line_c = shortest_line(line_1, line_2)

            bool_1 = line_1.boundary.touches(line_c)
            bool_2 = line_2.boundary.touches(line_c)

            if all([bool_1, bool_2]):
                line_1 = line_merge(MultiLineString([line_1, line_c]))

            else:
                if bool_1:
                    line_1 = line_merge(MultiLineString([line_1, line_c]))
                    line_2 = split(line_2, line_1)
                else:
                    line_1 = split(line_1, line_2)
                    line_2 = line_merge(MultiLineString([line_2, line_c]))

            line_dict[line_c.length] = {idx_1: line_1, idx_2: line_2}

        idx = min(line_dict.items())[1].keys()
        geoms = min(line_dict.items())[1].values()
        gdf_lines.loc[idx, "geometry"] = list(geoms)
        gdf_lines = gdf_lines.explode(index_parts=True)

    return gdf_lines


def graph_to_gdfs(graph):
    """
    Convert a nx.Graph to node and/or edge gpd.GeoDataFrame.
    inspired by osmnx package > utils_graph.graph_to_gdfs() function

    Parameters:
            graph (nx.Graph): graph representation of lines_gdf

    Returns:
            gpd.GeoDataFrames or tuple: gdf_nodes and gdf_edges as tuple
            of (gdf_nodes, gdf_edges). gdf_nodes is indexed by (x, y)
            coordinate and gdf_edges is multi-indexed by (u, v, key)
            following normal nx.Graph structure.
    """

    crs = graph.graph["crs"]

    if not graph.nodes:
        raise ValueError("Graph contains no nodes.")

    nodes, data = zip(*graph.nodes(data=True))
    geom = list(d["geometry"] for d in data)
    gdf_nodes = gpd.GeoDataFrame(data, index=nodes, crs=crs, geometry=geom)

    if not graph.edges:
        raise ValueError("Graph contains no edges.")

    u, v, data = zip(*graph.edges(data=True))
    geom = list(d["geometry"] for d in data)
    gdf_edges = gpd.GeoDataFrame(data, crs=crs, geometry=geom)

    gdf_edges["u"] = u
    gdf_edges["v"] = v
    gdf_edges.set_index(["u", "v"], inplace=True)

    return gdf_nodes, gdf_edges


def graph_from_gdfs(gdf_nodes, gdf_edges):
    """
    Convert node and edge GeoDataFrames to a Graph.
    This function is the inverse of `graph_to_gdfs` and is designed to work in
    conjunction with it.
    inspired by osmnx package > utils_graph.graph_from_gdfs() function

    Parameters
    ----------
    gdf_nodes : gpd.GeoDataFrame of graph nodes
    gdf_edges : gpd.GeoDataFrame of graph edges uniquely
                multi-indexed by u, v, key

    Returns
    -------
    graph : nx.Graph
    """

    graph_attrs = {"crs": gdf_edges.crs}
    graph = nx.Graph(**graph_attrs)

    gdf_edges.sort_index(inplace=True)
    attr_names = gdf_edges.columns.to_list()
    for (u, v), attr_vals in zip(gdf_edges.index, gdf_edges.to_numpy()):
        data_all = zip(attr_names, attr_vals)
        data = {
            name: val
            for name, val in data_all
            if isinstance(val, list) or pd.notna(val)
        }
        graph.add_edge(u, v, **data)

    gdf_nodes.sort_index(inplace=True)
    graph.add_nodes_from(set(gdf_nodes.index) - set(graph.nodes))
    for col in gdf_nodes.columns:
        nx.set_node_attributes(
            graph, name=col, values=gdf_nodes[col].dropna().to_dict()
        )

    return graph


def join_stations_to_graph(graph, gdf_stations):
    """
    Station coordinates are shifted to nearest point on rail lines. Lines
    are segmented based on nearest point coordinates. In case of duplicate
    nearest point allocation synthetic edges are added.

    Parameters:
            graph (nx.Graph): graph representation of lines_gdf

    Returns:
            gpd.GeoDataFrames or tuple: gdf_nodes and gdf_edges as tuple
            of (gdf_nodes, gdf_edges). Columns in graph gdfs contain all
            relevant attributes from lines_gdf / stations_gdf.
    """

    # prepare gdfs
    gdf_nodes, gdf_edges = graph_to_gdfs(graph)
    gdf_nodes["station"] = False
    gdf_edges["segmented"] = False
    gdf_edges["line_geom"] = gdf_edges.geometry

    # spatial join stations to rail lines,
    # remove stations outside search radius
    sjoin = gpd.sjoin_nearest(
        gdf_stations,
        gdf_edges,
        how="left",
        rsuffix="",
        distance_col="dist",
        max_distance=cfg["geo"]["search_radius"],
    ).sort_values("dist", ascending=True)
    sjoin = sjoin[sjoin.dist.notna()]
    sjoin = sjoin[~sjoin.index.duplicated(keep="first")]

    # for every station find nearest rail line coord,
    # union all points allocated to the same rail line
    sjoin["n_pt"] = sjoin.apply(
        lambda x: nearest_points(x.geometry, x.line_geom)[1], axis=1
    )
    sjoin["n_pts"] = sjoin.FID_rail_d.apply(
        lambda x: sjoin[sjoin.FID_rail_d == x].n_pt.unary_union
    )
    # split lines into segments based on n_pts (MultiPoint)
    sjoin["line_geom"] = sjoin.apply(
        lambda x: split(snap(x.line_geom, x.n_pts, 100), x.n_pts), axis=1
    )

    # rebuild gdf_nodes
    sjn_nodes = sjoin[gdf_stations.columns]
    sjn_nodes["station"] = True
    # station coords are shifted onto rail line
    sjn_nodes["geometry"] = sjoin["n_pt"]
    # when allocated duplicated, keep station coord and build connecting edges
    idx_dpl = sjn_nodes[sjn_nodes.geometry.duplicated()].index
    sjn_nodes.loc[idx_dpl, "geometry"] = gdf_stations.loc[idx_dpl, "geometry"]
    con_edges = sjoin.loc[idx_dpl, gdf_edges.columns.drop("line_geom")]
    con_edges["geometry"] = [
        LineString(geoms)
        for geoms in zip(
            gdf_stations.loc[idx_dpl, "geometry"], sjoin.loc[idx_dpl, "n_pt"]
        )
    ]
    # handle indexing
    # sjn_nodes = round_gdf(sjn_nodes)
    idx_nde = geom_to_int(sjn_nodes, "point")
    sjn_nodes = sjn_nodes.set_index(idx_nde)

    # rebuild gdf_edges
    sjn_edges = sjoin[gdf_edges.columns.drop("line_geom")]
    sjn_edges["segmented"] = True
    # overwrite gdf with segmented line geometries
    sjn_edges["geometry"] = sjoin["line_geom"]
    sjn_edges = sjn_edges.explode("geometry", index_parts=True)
    sjn_edges = sjn_edges.drop_duplicates("geometry", keep="first")
    # add connecting edges to frame
    sjn_edges = pd.concat([sjn_edges, con_edges], axis=0)
    # recalc geometries and length
    sjn_edges["geometry"] = sjn_edges.geometry
    sjn_edges["length"] = sjn_edges.geometry.length
    # handle indexing
    # sjn_edges = round_gdf(sjn_edges)
    sjn_edges.rename(columns={"index_0": "u", "index_1": "v"})
    sjn_edges["u"], sjn_edges["v"] = geom_to_int(sjn_edges, "line")
    sjn_edges = sjn_edges.set_index(["u", "v"])
    # remove original lines from gdf
    idx_drp = [tuple(idx) for idx in sjoin[["index_0", "index_1"]].values]
    gdf_edges = gdf_edges.sort_index().drop(idx_drp, axis=0)
    gdf_edges = gdf_edges.drop("line_geom", axis=1)

    # merge old and new gdfs
    gdf_nodes = pd.concat([gdf_nodes, sjn_nodes])
    gdf_edges = pd.concat([gdf_edges, sjn_edges])

    return gdf_nodes, gdf_edges


def remove_isolates(graph):
    largest_cc = max(nx.connected_components(graph), key=len)
    isolate_nodes = set(graph.nodes) - largest_cc
    graph.remove_nodes_from(isolate_nodes)

    # isolate_stations = [graph.nodes[n]["name"]
    # for n in isolate_nodes if graph.nodes[n]["station"]]
    logger.info(f"Removed {len(isolate_nodes)} isolated nodes from graph.")

    return graph


def remove_deadends(graph):
    while True:
        deadend_nodes = [
            n
            for n, deg in dict(graph.degree()).items()
            if deg == 1 and not graph.nodes[n]["station"]
        ]
        if len(deadend_nodes) > 0:
            graph.remove_nodes_from(deadend_nodes)
        else:
            break
    return graph
