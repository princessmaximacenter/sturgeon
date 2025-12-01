#!/usr/bin/env python3

import shutil
import time
import logging
import threading
import sys
from tokenize import String
import click
from functools import wraps
import signal
from pathlib import Path
import yaml
import json
from datetime import datetime
from watchdog.observers import Observer

from python_scripts import SturgeonBamHandling as SBH
from python_scripts import SturgeonLogging as SL
from python_scripts import SturgeonLivePlotting as SLP


# Load config file with default values
pythonPath = Path(__file__).resolve()
sturgeonRoot = pythonPath.parents[0]
CONFIG_PATH = sturgeonRoot / "config.yaml"

def load_config(config_yaml: Path) -> dict:
    """
    Load default config settings from yaml
    :param config_yaml: Path to config.yaml
    :return: dict containing default config settings
    """
    try:
        with open(config_yaml, 'r') as config_file:
            return yaml.safe_load(config_file)
    except Exception as e:
        logging.error(f"Error loading config: {e}, exiting...")
        sys.exit(1)

# Handle termination of python script
shutdown_event = threading.Event()

def _handle_exit(signum, frame) -> None:
    """
    Enters shutdown mode after receiving stop signal
    """
    SL._get_app_logger().info(f"Received termination signal ({signum}), shutting down...")
    shutdown_event.set()

def _register_signal_handlers():
    """
    Exiting shutdown either by killing through the CLI or via the GUI
    """
    signal.signal(signal.SIGINT, _handle_exit)
    signal.signal(signal.SIGTERM, _handle_exit)


def _normalize_barcode(ctx, param, value) -> String:
    """If users inputs 1 character string for barcode, add 0 as prefix"""
    if value and len(value) == 1:
        return f"0{value}"
    return value

class SturgeonPrediction:
    """
    Class to handle the entire process of a Sturgeon classification run
    -Setup
    -Classification and plotting
    -Cleanup and shutdown
    """

    def __init__(self, **kwargs):
        # Store all click parameters in attributes, also used for logging later
        for paramater, value in kwargs.items():
            setattr(self, paramater, value)

        self.app_log =  None
        self.lock_manager = None
        self.final_metadata = {}
        self.results_directory = None
        self.observer = None
        self.event_handler = None
        self.shutdown_event = shutdown_event
    def  _setup(self):
        """
        Handles output directory, logging initialization, and lock creation
        """

        # Verify output directory
        try:
            if not self.gui_activated:
                self.output.mkdir(parents=True, exist_ok=False)
            else:
                self.output.mkdir(parents=True, exist_ok=True)

        except FileExistsError:
            print(f"Error: Output directory '{self.output}' already exists. Exiting sturgeon...")
            sys.exit(1)

        # Initialization of logging
        SL._setup_logging(CONFIG_PATH, self.output)
        self.app_log = SL._get_app_logger()
        self.app_log.info(f"Output directory {self.output} created. Starting run")

        # Initialization of metadata output json
        all_cli_args = {k: getattr(self, k) for k in self.__dict__ if k not in ['app_log', 'lock_manager', 'final_metadata', 'results_directory', 'observer', 'event_handler', 'shutdown_event']}
        loggable_params = {
            k: str(v.resolve()) if isinstance(v, Path) else v
            for k, v in all_cli_args.items()
        }
        self.final_metadata = {
            "RunParamaters": [loggable_params],
            "RunInfo": [
                {
                    "run_time_start": datetime.now().isoformat()
                }
            ],
            "RunResults": []
        }

        # Finalize setup with handling of results directory and lock file
        self.results_directory = self._set_results_directory()
        _register_signal_handlers()
        self.lock_manager = SBH.LockManager(self.lock)
        self.lock_manager._check_lock()
        self.app_log.info("Setup complete. Ready for execution")

    def _set_results_directory(self) -> Path:
        """
        Finalize the output directory for sturgeon results
        :return: Path to directory with bam files to watched and processed
        """
        if self.gridion:
            return self.input
        else:
            if self.barcode != "unclassified":
                self.barcode = f"barcode{self.barcode}"
            else:
                pass
            base_input = self.input.resolve()
            # Option 1: ends with bam_pass/barcode → use as-is
            if base_input.parts[-2:] == ("bam_pass", self.barcode):
                return base_input

            # Option 2: ends with bam_pass → add barcode
            elif base_input.name == "bam_pass":
                return Path(f"{base_input}/{self.barcode}")

            # Option 3: assume it's the sequencing dir → add bam_pass/barcode
            else:
                return Path(f"{base_input}/bam_pass/{self.barcode}")


    def _wait_for_input_directory(self) -> bool:
        """
        Waits for the results directory with the BAM files to be created.
        Checks every 1 second, but only prints every 30 to keep the log clean.
        Returns False if shutdown is requested.
        """
        wait_time = 0
        log_interval_time = 30

        while not self.shutdown_event.is_set():
            if self.shutdown_file.exists():
                self.app_log.info("Shutdown requested during wait for input directory.")
                self.shutdown_event.set()
                return False
            if self.input.exists():
                self.app_log.info(f"Results directory {self.input} found, proceeding with sturgeon analysis.")
                return True
            if wait_time % log_interval_time == 0:
                self.app_log.info(f"Waiting for results directory {self.input} to be created.")
            time.sleep(1)
            wait_time += 1

        self.app_log.info("Shutdown event set before input directory appeared.")
        return False

    def _file_cleanup_shutdown(self) -> None:
        """
        After shutdown of sturgeon move last created analysis files to main output folder.
        Also cleans up symbolic links in merged_bam dir to avoid potential problems later.
        """

        last_iteration = self.event_handler.iteration
        final_dir = Path(f"{self.output}/iteration_{last_iteration}")

        if not final_dir.is_dir():
            last_iteration -= 1
            final_dir = Path(f"{self.output}/iteration_{last_iteration}")
            if not Path(f"{final_dir}/CNV_plot_iteration_{last_iteration}.pdf").is_file():
                self.event_handler.plot_cnv(last_iteration)
        else:
            if not Path(f"{final_dir}/CNV_plot_iteration_{last_iteration}.pdf").is_file():
                self.event_handler.plot_cnv()
            if not Path(f"{final_dir}/classifier_progress_iteration_{last_iteration}.tsv").is_file():
                self.event_handler.plot_process()

        try:
            src_files = [f"CNV_plot_iteration_{last_iteration}.pdf", f"confidence_over_time_plot_iteration_{last_iteration}.pdf",
                         f"merged_probes_methyl_calls_*_iteration_{last_iteration}.pdf",
                         f"classifier_progress_iteration_{last_iteration}.tsv"]
            final_classification_dest = Path(f"{self.output}/")
            for pattern in src_files:
                for final_classification_src in Path(final_dir).glob(pattern):
                    shutil.copy(final_classification_src, final_classification_dest)

        except Exception as e:
            self.app_log.warning(f"Could not copy final classification file: {e}")

        #Remove symbolic links to avoid potential problems in the future
        bam_link_dir = Path(f"{self.output}/merged_bams/")
        for bam_link in bam_link_dir.glob("bam_for_CNV_it*"):
            bam_link.unlink()

    def _cleanup(self):
        """
        Handles shutdown of wrapper script, final file handling and lock/shutdown flag removal
        """

        self.app_log.info("Starting cleanup process...")
        if self.observer:
            self.observer.stop()
            self.observer.join()

        if self.lock_manager:
            self.lock_manager._remove_lock()

        if self.event_handler:
            self.app_log.info("Moving final plots and cleaning up files before shutting down Sturgeon")
            self._file_cleanup_shutdown()

            try:
                classification_data = SLP.get_final_classification(self.output, self.event_handler.iteration)
                self.final_metadata["RunResults"] = [classification_data]
                self.app_log.info(f"Final results logged: Class={classification_data['final_classification']}, Score={classification_data['final_score']}")

            except Exception as e:
                self.app_log.error(f"Failed to log final classification metadata: {e}", exc_info=True)
                self.final_metadata["RunResults"] = [{"error": str(e)}]

        metadata_file = f"{self.output}/sturgeon_metadata.json"
        try:
            with open(metadata_file, 'w') as f:
                json.dump(self.final_metadata, f, indent=2)
            self.app_log.info(f"Final metadata report saved to {metadata_file}")
        except Exception as e:
            self.app_log.error(f"Failed to write final metadata JSON: {e}")

        if self.shutdown_file.exists():
            self.shutdown_file.unlink()
            self.app_log.info("Removed shutdown flag")
        else:
            self.app_log.info("No shutdown flag found to remove")

    def run(self):
        """
        Main execution loop
        Calls setup, then execution, then cleanup
        """

        # Setup
        try:
            self._setup()
        except SystemExit:
            # Clean exit if setup for some reason failed
            return
        except Exception as e:
            print(f"Fatal ERROR during setup: {e}")
            sys.exit(1)

        # Execution main loop
        try:
            if not self.results_directory.exists():
                if not self._wait_for_input_directory():
                    return

            self.app_log.info(f"Starting to monitor for new BAM files in: {self.results_directory}")
            self.event_handler = SBH.NewBamFileHandler(
                self.sturgeon_script, self.output, self.model, self.freq, self.utils, self.r_script,
                self.results_directory, self.gridion, self.shutdown_event, self.live_run, self.version2, self.conf)
            self.observer = Observer()
            self.observer.schedule(self.event_handler, path=self.results_directory, recursive=False)
            self.observer.start()

            while not self.shutdown_event.is_set():
                if self.shutdown_file.exists():
                    self.app_log.info("Shutdown file detected. Initiating shutdown...")
                    self.app_log.info("Waiting for running process to complete")
                    self.event_handler.wait_for_process_completion()
                    self.shutdown_event.set()
                time.sleep(1)

            self.app_log.info("Shutdown requested. Cleaning up...")

            # Final Cleanup
        finally:
            self._cleanup()


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
        "-b", "--barcode", type=str, default=None, callback=_normalize_barcode,help="Barcode used in library preparation."
    )
    @click.option(
        "-f", "--freq", type=int, default=None, help="Number of iterations before merging BAMs and plotting CNV."
    )
    @click.option(
        "-m", "--model", type=click.Path(path_type=Path, exists=True), default=None, help="Location of model used for sturgeon prediction"
    )
    @click.option(
        "-u", "--utils", type=click.Path(path_type=Path, exists=True), default=None, help="Location of utils directory"
    )
    @click.option(
        "-r", "--r_script", type=click.Path(path_type=Path, exists=True), default=None,help="Location of R script for plotting CNV"
    )
    @click.option(
        "-g", "--gridion", type=bool, default = False, help = "If run is a gridion verification run, some parameters are changed"
    )
    @click.option(
        "-sf", "--shutdown_file", type=click.Path(path_type=Path, exists=False), default=None,help="Location of shutdown flag"
    )
    @click.option(
        "-v2", "--version2", type=str, default=False, help = "If true, prediction will be run with sturgeon model v2"
    )
    @click.option(
        "--gui_activated", is_flag=True, default=False, help="Flag to indicate script is run through GUI"
    )
    @click.option(
        "-lr", "--live_run", is_flag=True, default=False, help="Flag to indicate whether sequencing and analysis is live"
    )
    @click.option(
        "--conf", type=click.Path(path_type=Path, exists=True), default=None, help="Location of python script for confidence over time for sturgeon model V2"
    )

    @wraps(func)
    def wrapper(*args, **kwargs):
        config = load_config(CONFIG_PATH)
        kwargs['input'] = kwargs.get('input') or Path(config['paths']['bam_input'])
        kwargs['output'] = kwargs.get('output') or Path(config['paths']['results_dir'])
        kwargs['lock'] = kwargs.get('lock') or Path(config['paths']['script_lock'])
        kwargs['sturgeon_script'] = kwargs.get('sturgeon_script') or Path(config['paths']['sturgeon_script'])
        kwargs['model'] = kwargs.get('model') or Path(config['paths']['model'])
        kwargs['utils'] = kwargs.get('utils') or Path(config['paths']['utils'])
        kwargs['r_script'] = kwargs.get('r_script') or Path(config['paths']['r_script'])
        kwargs['barcode'] = kwargs['barcode'] or config.get('barcode')
        kwargs["freq"] = kwargs["freq"] or config.get("freq")
        kwargs["gridion"] = kwargs["gridion"] or config.get("gridion")
        kwargs["shutdown_file"] = kwargs["shutdown_file"] or Path(config['paths']['shutdown_file'])
        kwargs["version2"] = kwargs["version2"] or config.get("sturgeon_V2")
        kwargs["gui_activated"] = kwargs["gui_activated"] or config.get("gui_activated")
        kwargs["live_run"] = kwargs["live_run"] or config.get("live_run")
        kwargs["conf"] = kwargs["conf"] or Path(config['paths']["conf_time_plot"])

        return func(*args, **kwargs)

    return wrapper


@click_command
def main(**kwargs):
    """
    Sturgeon Prediction:
    Monitors a directory for new BAM files (live run) or processes all BAM files (post-sequencing run) and performs Sturgeon classification
    """

    runner = SturgeonPrediction(**kwargs)
    runner.run()

if __name__ == "__main__":
    main()