import logging

import geopandas as gpd
import yaml
from shapely import line_merge, wkt
from shapely.geometry import GeometryCollection, MultiLineString
from shapely.ops import snap, split


def read_config():
    with open("config.yaml") as ymlfile:
        cfg = yaml.safe_load(ymlfile)

    return cfg


cfg = read_config()


def setup_logger(name):
    formatter = logging.Formatter(cfg["logger"]["format"])

    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger


def round_gdf(gdf):
    precision = cfg["geo"]["precision"]

    geo_cols = gdf.select_dtypes(include="geometry").columns
    for col in geo_cols:
        geoms = gdf.apply(lambda x: wkt.dumps(x[col], True, precision), axis=1)
        gdf[col] = gpd.GeoSeries.from_wkt(geoms)

    return gdf


def tuple_coords(coords):
    return tuple(int(c) for c in coords)


def geom_to_int(gdf, geom_type, col="geometry"):
    if geom_type == "point":
        s = gdf.apply(lambda x: tuple_coords((x[col].x, x[col].y)), axis=1)
        return s

    elif geom_type == "line":
        s1 = gdf.apply(lambda x: tuple_coords(x[col].coords[0]), axis=1)
        s2 = gdf.apply(lambda x: tuple_coords(x[col].coords[-1]), axis=1)
        return [s1, s2]

    else:
        TypeError("Please refer to the format for parameter 'geom_type'.")


def assemble_line(line_12, line_c, deg_12):
    if deg_12 == 1:
        line_12 = line_merge(MultiLineString([line_12, line_c]))
    else:
        if line_12.geom_type == "LineString":
            line_12 = GeometryCollection([line_12, line_c])
        else:
            geoms = list(line_12.geoms) + [line_c]
            line_12 = GeometryCollection(geoms)

    return line_12


def split_line(line, splitter):
    tol = cfg["geo"]["precision"]
    line_s = split(snap(line, splitter, tol), splitter)

    return line_s
