# SturgeonLogging.py
import logging
import logging.config
import yaml
import json
from pathlib import Path


def _get_app_logger():
    return logging.getLogger('sturgeon.app')

def _setup_logging(config_path: Path, output_dir: Path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
        log_config = config.get('logging')
        if not log_config:
            raise ValueError("YAML config must contain a logging section")

        app_log_path = output_dir / 'sturgeon_run.log'
        if 'file_app' in log_config['handlers']:
            # Replace the placeholder in the config with the final path string
            log_config['handlers']['file_app']['filename'] = str(app_log_path)

        logging.config.dictConfig(log_config)