!!! Active development - not usable at the moment !!!
# railwaynetworks

A small package to build simplified graph representations of existing railway 
track networks. Scope of this tool is to avoid the pre-processing of large 
[OSM data](https://download.geofabrik.de/) by using the following open data sources instead
- [DIVA GIS](https://www.diva-gis.org/)
- [trainline-eu/stations](https://github.com/trainline-eu/stations)

The code has been developed and tested for the European continent, but is 
theoretically applicable / adaptable for the whole world.

Please note that created network lines in general do not contain operational information
(track usages, gauges or max speeds) and some included railway tracks may be out of service.

## Installation
create, activate & install virtual environment
```commandline
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
