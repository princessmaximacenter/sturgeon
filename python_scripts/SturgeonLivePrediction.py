import shutil
import time
import logging
import threading
import sys
from tokenize import String
import click
import signal
from pathlib import Path
import yaml
from watchdog.observers import Observer

from python_scripts import SturgeonBamHandling as SBH
from python_scripts import SturgeonLogging as SL


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

SL._setup_logging(CONFIG_PATH)
log = SL._get_logger()


# Handle termination of python script
shutdown_event = threading.Event()

def _handle_exit(signum, frame) -> None:
    """
    Enters shutdown mode after receiving stop signal
    """
    log.info(f"Received termination signal ({signum}), shutting down...")
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

def set_results_directory(input: Path, barcode: str, gridion: bool) -> Path:
    """
    Finalize the output directory for sturgeon results
    :param input: Results path of sequencing run
    :param barcode: Given barcode used for library prep
    :param gridion: Boolean for checking if run is a gridion verification run
    :return: Path to directory with bam files to watched and processed
    """
    if gridion:
        return input
    else:
        if barcode != "unclassified":
            barcode = f"barcode{barcode}"
        else:
            pass
        base_input = input.resolve()
        # Option 1: ends with bam_pass/barcode → use as-is
        if base_input.parts[-2:] == ("bam_pass", barcode):
            return base_input

        # Option 2: ends with bam_pass → add barcode
        elif base_input.name == "bam_pass":
            return Path(f"{base_input}/{barcode}")

        # Option 3: assume it's the sequencing dir → add bam_pass/barcode
        else:
            return Path(f"{base_input}/bam_pass/{barcode}")


def _wait_for_input_directory(input: Path, shutdown_file: Path) -> bool:
    """
    Waits for the results directory with the BAM files to be created.
    Checks every 1 second, but only prints every 30 to keep the log clean.
    Returns False if shutdown is requested.
    """
    wait_time = 0
    log_interval_time = 30

    while not shutdown_event.is_set():
        if shutdown_file.exists():
            log.info("Shutdown requested during wait for input directory.")
            shutdown_event.set()
            return False
        if input.exists():
            log.info(f"Results directory {input} found, proceeding with sturgeon analysis.")
            return True
        if wait_time % log_interval_time == 0:
            log.info(f"Waiting for results directory {input} to be created.")
        time.sleep(1)
        wait_time += 1

    log.info("Shutdown event set before input directory appeared.")
    return False

def _file_cleanup_shutdown(output, iteration, plot_cnv, plot_process) -> None:
    """
    After shutdown of sturgeon move last created analysis files to main output folder.
    Also cleans up symbolic links in merged_bam dir to avoid potential problems later.

    :param output: Results directory
    :param iteration: Last iteration of analysis
    :param plot_cnv: function to create cnv plot
    :param plot_process: function to create the tsv file with classifier scores per iteration
    """

    final_dir = Path(f"{output}/iteration_{iteration}")
    if not final_dir.is_dir():
        iteration -= 1
        final_dir = Path(f"{output}/iteration_{iteration}")
        if not Path(f"{final_dir}/CNV_plot_iteration_{iteration}.pdf").is_file():
            plot_cnv(iteration)
    else:
        if not Path(f"{final_dir}/CNV_plot_iteration_{iteration}.pdf").is_file():
            plot_cnv()
        if not Path(f"{final_dir}/classifier_progress_iteration_{iteration}.tsv").is_file():
            plot_process()

    try:
        src_files = [f"CNV_plot_iteration_{iteration}.pdf", f"confidence_over_time_plot_iteration_{iteration}.pdf",
                     f"merged_probes_methyl_calls_*_iteration_{iteration}.pdf",
                     f"classifier_progress_iteration_{iteration}.tsv"]
        final_classification_dest = Path(f"{output}/")
        for pattern in src_files:
            for file_name in Path(final_dir).glob(pattern):
                shutil.copy(file_name, final_classification_dest)

    except Exception as e:
        log.warning(f"Could not copy final classification file: {e}")

    #Remove symbolic links to avoid potential problems in the future
    bam_link_dir = Path(f"{output}/merged_bams/")
    for bam_link in bam_link_dir.glob("bam_for_CNV_it*"):
        bam_link.unlink()


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
        "-v2", "--version2", type=str, default = False, help = "If true, prediction will be run with sturgeon model v2"
    )
    @click.option(
        "--gui_activated", is_flag=True, default=False, help="Flag to indicate script is run through GUI"
    )
    @click.option(
        "--conf", type=click.Path(path_type=Path, exists=True), default=None, help="Location of python script for confidence over time for sturgeon model V2"
    )

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
        kwargs["conf"] = kwargs["conf"] or Path(config['paths']["conf_time_plot"])
        return func(*args, **kwargs)

    return wrapper


@click_command
def main(input: Path, output: Path, lock: Path, sturgeon_script: Path, barcode: str, freq: int, model: Path, utils: Path, r_script: Path, gridion: bool, shutdown_file:Path, version2: str, gui_activated: bool, conf: Path):
    """
    Monitors a folder for new BAM files and processes them as they appear.
    """

    results_directory = set_results_directory(input, barcode, gridion)
    _register_signal_handlers()

    lock_manager = SBH.LockManager(lock)
    lock_manager._check_lock()

    observer = None
    event_handler = None

    try:
        if not gui_activated:
            output.mkdir(exist_ok=False)  # Check if output directory exists
        else:
            log.info("Sturgeon started through gui, output dir already created")
    except FileExistsError:
        log.warning("Output directory already exists, exiting sturgeon...")

    else:
        if not results_directory.exists():
            if not _wait_for_input_directory(results_directory, shutdown_file):
                return

        log.info(f"Starting to monitor for new BAM files in: {results_directory}")
        event_handler = SBH.NewBamFileHandler(sturgeon_script, output, model, freq, utils, r_script, results_directory,
                                              gridion, shutdown_event, version2, conf)
        observer = Observer()
        observer.schedule(event_handler, path=results_directory, recursive=False)
        observer.start()

        while not shutdown_event.is_set():
            if shutdown_file.exists():
                log.info("Shutdown file detected. Initiating shutdown...")
                event_handler.wait_for_process_completion()
                shutdown_event.set()
            time.sleep(1)

        log.info("Shutdown requested. Cleaning up...")
    finally:
        if observer:
            observer.stop()
            observer.join()

        lock_manager._remove_lock()

        if event_handler:
            log.info("Moving final plots and cleaning up files before shutting down Sturgeon")
            _file_cleanup_shutdown(output, event_handler.iteration, event_handler.plot_cnv, event_handler.plot_process)

        if shutdown_file.exists():
            shutdown_file.unlink()
            log.info("Removed shutdown flag")
        else:
            log.info("No shutdown flag found to remove")


if __name__ == "__main__":
    main()