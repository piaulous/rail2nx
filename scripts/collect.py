import io
import logging
import os.path
import zipfile

import geopandas as gpd
import pandas as pd
import requests
from shapely.geometry import Point

from scripts.tools import read_config

cfg = read_config()
logger = logging.getLogger("rail2nx")


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

    shps_dir = cfg["dirs"]["lines"]
    cols_raw = ["FID_rail_d", "EXS_DESCRI", "FCO_DESCRI", "ISO", "geometry"]
    cols_tar = ["rail_id", "status", "tracks", "country", "geometry"]
    ccodes_a2 = get_country_codes(scope, "alpha-2")
    ccodes_a3 = get_country_codes(scope, "alpha-3")

    logger.info("Collect rail line data from DIVA-GIS")

    rlines = []
    for cc in ccodes_a3:
        cc_fname = f"{shps_dir}/{cc}_rails.shp"

        if not os.path.isfile(cc_fname):
            diva_url = f"https://biogeo.ucdavis.edu/data/diva/rrd/{cc}_rrd.zip"
            r = requests.get(diva_url)
            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(shps_dir)

        # bugfix workaround
        if cc == "ROU":
            cc_fname = f"{shps_dir}/ROM_rails.shp"

        cc_rlines = gpd.read_file(cc_fname, crs=cfg["geo"]["crs_def"]).to_crs(
            cfg["geo"]["crs_eur"]
        )[cols_raw]
        rlines.append(cc_rlines)

    cc_mapping = dict(zip(ccodes_a3, ccodes_a2))
    cc_mapping.update({"ROM": "RO"})  # bugfix workaround
    col_mapping = dict(zip(cols_raw, cols_tar))

    gdf_lines = pd.concat(rlines).reset_index(drop=True)
    gdf_lines = gdf_lines.rename(columns=col_mapping)
    gdf_lines["country"] = gdf_lines.country.replace(cc_mapping, regex=True)

    return gdf_lines


def get_rail_stations(scope, locate_coords=True):
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

    fname = cfg["dirs"]["stations"]

    if not os.path.isfile(fname):
        stations = pd.read_csv(
            cfg["urls"]["stations"],
            delimiter=";",
            low_memory=False,
        )
        stations.to_csv(fname, index=False)

    stations = pd.read_csv(fname, low_memory=False)
    x, y = stations.longitude, stations.latitude
    stations["geometry"] = gpd.points_from_xy(x, y)
    stations = (
        gpd.GeoDataFrame(stations, crs=cfg["geo"]["crs_def"])
        .set_geometry("geometry")
        .to_crs(cfg["geo"]["crs_eur"])
    )

    # remove stations with duplicated coordinates
    stations = stations.sort_values("parent_station_id", ascending=True)
    stations = stations[
        ~stations.duplicated(subset=["latitude", "longitude"], keep="last")
    ]

    rail_operators = [
        "sncf_id",
        "entur_id",
        "db_id",
        "cff_id",
        "leoexpress_id",
        "obb_id",
        "ouigo_id",
        "trenitalia_id",
        "trenord_id",
        "ntv_id",
        "hkx_id",
        "renfe_id",
        "atoc_id",
        "benerail_id",
        "westbahn_id",
    ]

    stations["has_rail_id"] = stations[rail_operators].notnull().any(axis=1)
    stations = stations.loc[stations["has_rail_id"]]
    bus_filter = "bussteig|bussterminal|autobus"
    stations = stations[~stations["slug"].str.contains(bus_filter, na=False)]

    if isinstance(scope, str):
        if len(scope) == 2:
            stations = stations[stations.country == scope]
        else:
            stations = stations[stations.time_zone.str.contains(scope)]
    elif isinstance(scope, list):
        stations = stations[stations.country.isin(scope)]
    else:
        TypeError("Please refer to the format for parameter 'scope'.")

    if locate_coords:
        locate_missing_coords(stations)

    return stations


def locate_missing_coords(stations_gdf):
    fname = cfg["dirs"]["stations"]
    nom = cfg["nominatim"]
    stations = stations_gdf[stations_gdf.geometry.is_empty].set_crs(
        cfg["geo"]["crs_def"], allow_override=True
    )

    unlocatables = []

    for idx, row in stations.iterrows():
        logger.info(f"Trying to locate station '{row['name']}'")

        query = (
            f"{nom['url']}&"
            f"accept-language={nom['language']}&"
            f"format={nom['format']}&"
            f"q={row.slug}+{nom['tags']}&"
            f"countrycodes={row.country}&"
            f"addressdetails=1"
        )
        response = requests.get(query)

        if response.json():
            is_city = "f"  # found station
            osm = response.json()[0]

        else:
            query = query.replace(
                f"+{nom['tags']}", ""
            )  # remove railway / station tag from query
            response = requests.get(query)
            matched_entries = [
                entry for entry in response.json() if entry["class"] == "place"
            ]

            if matched_entries:
                is_city = "t"  # place with station's name has been found
                osm = matched_entries[0]

            else:
                logger.info(
                    f"Could not locate station '{row['name']}', "
                    f"entry will be removed."
                )
                unlocatables.append(idx)
                continue

            logger.info(f"Successfully located station '{row['name']}'")

        # fill missing data
        stations.at[idx, "latitude"] = float(osm["lat"])
        stations.at[idx, "longitude"] = float(osm["lon"])
        stations.at[idx, "geometry"] = Point(osm["lon"], osm["lat"])
        stations.at[idx, "is_city"] = is_city

    stations = stations.to_crs(cfg["geo"]["crs_eur"])
    stations_gdf.update(stations)

    # update coords and drop stations that could not be located
    if len(stations):
        stations_csv = pd.read_csv(fname, low_memory=False)
        stations_csv.update(stations)
        stations_csv.drop(unlocatables, axis=0, inplace=True)
        stations_csv.to_csv(fname, index=False)
