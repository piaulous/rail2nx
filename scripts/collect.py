import io
import logging
import os.path
import zipfile

import geopandas as gpd
import pandas as pd
import requests

from scripts.tools import read_config

cfg = read_config()
logger = logging.getLogger("railwaynetworks")


def get_country_codes(scope, format="alpha-3"):
    """
    Returns country codes for specified scope in 2-letter or 3-letter format.

    Parameters:
            scope (str or list):
                continents ('Europe'),
                single country ('DE') or
                country list (['BE', 'NL', ...])
            format (str): 'alpha-2' or 'alpha-3'

    Returns:
            ccodes (pd.Series): country codes
    """

    fname = "data/country_codes.csv"

    if not os.path.isfile(fname):
        ccodes_url = cfg["urls"]["ccodes"]
        ccodes = pd.read_csv(ccodes_url)
        ccodes.to_csv(fname)

    ccodes = pd.read_csv(fname)

    if isinstance(scope, str):
        if len(scope) == 2:
            ccodes = ccodes[ccodes["alpha-2"] == scope]
        else:
            ccodes = ccodes[ccodes.region == scope]
    elif isinstance(scope, list):
        ccodes = ccodes[ccodes["alpha-2"].isin(scope)]
    else:
        TypeError("Please refer to the format for parameter 'scope'.")

    logger.info(f"Following countries are included {list(ccodes['alpha-2'])}")

    return ccodes[format]


def get_rail_lines(scope):
    """
    Country-wise download or import of railroads from DIVA-GIS database
    for specified scope.

    Source:
            https://www.diva-gis.org/gdata

    Parameters:
            scope (str or list):
                continents ('Europe')
                single country ('DE') or
                country list (['BE', 'NL', ...])

    Returns:
            lines_gdf (gpd.GeoDataFrame): LineStrings plus additional metadata
    """

    shps_dir = "data/lines_raw"
    ccodes = get_country_codes(scope)
    rlines = []

    logger.info("Collect rail line data from DIVA-GIS")

    for cc in ccodes:
        cc_fname = f"{shps_dir}/{cc}_rails.shp"

        if not os.path.isfile(cc_fname):
            diva_url = f"https://biogeo.ucdavis.edu/data/diva/rrd/{cc}_rrd.zip"
            r = requests.get(diva_url)
            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(shps_dir)

        if cc == "ROU":
            cc_fname = f"{shps_dir}/ROM_rails.shp"

        cc_rlines = gpd.read_file(cc_fname, crs=cfg["proj"]["crs_def"]).to_crs(
            cfg["proj"]["crs_eur"]
        )
        rlines.append(cc_rlines)

    return pd.concat(rlines)


def get_rail_stations(scope):
    """
    Download or import of trainline-eu/stations database for specified scope.

    Source:
            https://github.com/trainline-eu/stations

    Parameters:
            scope (str or list):
                continents ('Europe'),
                single country ('DE') or
                country list (['BE', 'NL', ...])

    Returns:
            stations_gdf (gpd.GeoDataFrame): stations with associated metadata
    """

    logger.info("Collect station data from trainline-eu/stations")

    fname = "data/stations_raw.csv"

    if not os.path.isfile(fname):
        stations = pd.read_csv(
            cfg["urls"]["stations"],
            delimiter=";",
            low_memory=False,
        )
        stations.to_csv(fname)

    stations = pd.read_csv(fname, low_memory=False)
    x, y = stations.longitude, stations.latitude
    stations["geometry"] = gpd.points_from_xy(x, y)
    stations = (
        gpd.GeoDataFrame(stations, crs=cfg["proj"]["crs_def"])
        .set_geometry("geometry")
        .to_crs(cfg["proj"]["crs_eur"])
    )

    if isinstance(scope, str):
        if len(scope) == 2:
            stations = stations[stations.country == scope]
        else:
            stations = stations[stations.time_zone.str.contains(scope)]
    elif isinstance(scope, list):
        stations = stations[stations.country.isin(scope)]
    else:
        TypeError("Please refer to the format for parameter 'scope'.")

    return stations
