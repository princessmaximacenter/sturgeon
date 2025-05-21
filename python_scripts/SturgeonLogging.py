# SturgeonLogging.py
import logging
import yaml
import sys
from pathlib import Path

def _get_logger():
    return logging.getLogger('root')

def _setup_logging(config_path: Path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
        log_level = config['logging']['log_level'].upper()

        logging.basicConfig(
            format='| %(levelname)-8s | %(message)s',
            level=logging.getLevelName(log_level),
            handlers=[
                logging.StreamHandler(sys.stdout)  # Only stream to stdout
            ],
            force=True
        )

