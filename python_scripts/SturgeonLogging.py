# SturgeonLogging.py
import logging
import logging.config
import yaml
import json
from pathlib import Path

class MetadataJsonFormatter(logging.Formatter):
    """
    A custom formatter to transform log record into a json file
    """
    def format(self, record):
        log_record = {
            "timestamp": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "logger": record.name,
            "module": record.module,
            "funcName": record.funcName,
            "lineno": record.lineno,
        }
        # Check if the message is the structured data we want to log
        if isinstance(record.msg, dict):
            # Merge the dictionary message into the log record
            log_record.update(record.msg)
            # Remove the 'message' key if it exists, as it's redundant
            if 'message' in log_record:
                del log_record['message']
        else:
            # Fallback for plain string messages logged to this handler
            log_record["message"] = record.getMessage()

        return json.dumps(log_record)


def _get_app_logger():
    return logging.getLogger('sturgeon.app')

def _get_metadata_logger():
    return logging.getLogger('sturgeon.metadata')




def _setup_logging(config_path: Path, output_dir: Path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

        log_config = config.get('logging')
        if not log_config:
            raise ValueError("YAML config must contain a logging section")

        metadata_filepath = str(output_dir / 'sturgeon_metadata.json')
        placeholder = "PLACEHOLDER STRING"
        handlers_config = log_config.get('handlers', {})

        if 'metadata' in handlers_config:
            handler = handlers_config['metadata']
            if handler.get('filename') == placeholder:
                handler['filename'] = metadata_filepath
            else:
                logging.getLogger('root').warning(f"Logging placeholder '{placeholder}' not found in config")


        logging.config.dictConfig(log_config)