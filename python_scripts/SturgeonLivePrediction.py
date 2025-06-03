import time
import logging
import threading
import sys
import click
import signal
from pathlib import Path
import yaml
from watchdog.observers import Observer

import SturgeonBamHandling as SBH
import SturgeonLogging as SL

# Load config file with default values
CONFIG_PATH = Path("/Users/k.v.cammel/Documents/Projects/P2_wrapper/sturgeon/python_scripts/config.yaml")
def load_config(config_yaml: Path) -> dict:
    """Load configuration from YAML"""
    try:
        with open(config_yaml, 'r') as config_file:
            return yaml.safe_load(config_file)
    except Exception as e:
        logging.error(f"Error loading config: {e}, exiting...")
        sys.exit(1)

SL._setup_logging(CONFIG_PATH)
log = SL._get_logger()


# Handle termination of python script
shutdown_event = threading.Event()

def handle_exit(signum, frame):
    log.info(f"Received termination signal ({signum}), shutting down...")
    shutdown_event.set()

def register_signal_handlers():
    signal.signal(signal.SIGINT, handle_exit)
    signal.signal(signal.SIGTERM, handle_exit)

def wait_for_input_directory(input: Path) -> None:
    """Wait for the results directory to be created"""
    while not input.exists():
        log.info(f"Waiting for results directory {input} to be created.")
        time.sleep(5)
        if shutdown_event.is_set():
            return
    log.info(f"Results directory {input} found, proceeding with sturgeon analysis.")

# Custom decorator for click command
def click_command(func):
    @click.command()
    @click.option(
        "-i", "--input", type=click.Path(path_type=Path, exists=False, dir_okay=True), default=None, help="Directory of sequencing run for sturgeon analysis"
    )
    @click.option(
        "-o", "--output", type=click.Path(path_type=Path, exists=False, writable=True), default=None, help="Directory where results are written."
    )
    @click.option(
        "-l", "--lock", type=click.Path(path_type=Path, exists=False, writable=True), default=None, help="Name of lock file."
    )
    @click.option(
        "-s", "--sturgeon_script", type=click.Path(path_type=Path, exists=True, file_okay=True), default=None, help="Path to the script that will be called for processing."
    )
    @click.option(
        "-b", "--barcode", type=int, default=None, help="Barcode used in library preparation."
    )
    @click.option(
        "-f", "--freq", type=int, default=None, help="Number of iterations before merging BAMs and plotting CNV."
    )
    @click.option(
        "-m", "--model", type=click.Path(path_type=Path, exists=True), default=None, help="Location of model used for sturgeon prediction"
    )
    @click.option(
        "-r", "--r_script", type=click.Path(path_type=Path, exists=True), default=None,help="Location of R script for plotting CNV"

    )
    @click.option(
        "-sf", "--shutdown_file", type=click.Path(path_type=Path, exists=False), default=None,help="Location of shutdown flag"
    )

    def wrapper(*args, **kwargs):
        config = load_config(CONFIG_PATH)
        kwargs['input'] = kwargs.get('input') or Path(config['paths']['bam_input'])
        kwargs['output'] = kwargs.get('output') or Path(config['paths']['results_dir'])
        kwargs['lock'] = kwargs.get('lock') or Path(config['paths']['script_lock'])
        kwargs['sturgeon_script'] = kwargs.get('sturgeon_script') or Path(config['paths']['sturgeon_script'])
        kwargs['model'] = kwargs.get('model') or Path(config['paths']['model'])
        kwargs['r_script'] = kwargs.get('r_script') or Path(config['paths']['r_script'])
        kwargs['barcode'] = kwargs['barcode'] or config.get('barcode')
        kwargs["freq"] = kwargs["freq"] or config.get("freq")
        kwargs["shutdown_file"] = kwargs["shutdown_file"] or Path(config['paths']['shutdown_file'])

        return func(*args, **kwargs)

    return wrapper


@click_command
def main(input: Path, output: Path, lock: Path, sturgeon_script: Path, barcode: int, freq: int, model: Path, r_script: Path, shutdown_file:Path):
    """
    Monitors a folder for new BAM files and processes them as they appear.
    """

    results_directory = Path(f"{input}/bam_pass/barcode{barcode}/")

    register_signal_handlers()

    lock_manager = SBH.LockManager(lock)
    lock_manager.check_lock()

    try:

        if not results_directory.exists():
            wait_for_input_directory(results_directory)
        if shutdown_event.is_set():
            log.info("Sturgeon terminated during wait for results")
            return
        log.info(f"Starting to monitor for new BAM files in: {results_directory}")
        event_handler = SBH.NewBamFileHandler(sturgeon_script, output, model, freq, r_script, results_directory)
        observer = Observer()
        observer.schedule(event_handler, path=results_directory, recursive=False)
        observer.start()

        while not shutdown_event.is_set():
            if shutdown_file.exists():
                log.info("Shutdown file detected. Initiating shutdown...")
                shutdown_event.set()
            time.sleep(1)

        log.info("Shutdown requested. Cleaning up...")
    finally:
        observer.stop()
        observer.join()
        lock_manager.remove_lock()

        if shutdown_file.exists():
            shutdown_file.unlink()
            log.info("Removed shutdown flag")
if __name__ == "__main__":
    main()