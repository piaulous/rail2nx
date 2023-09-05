import logging

import yaml


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
