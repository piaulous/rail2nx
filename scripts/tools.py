import logging

import geopandas as gpd
import yaml
from shapely import wkt


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
