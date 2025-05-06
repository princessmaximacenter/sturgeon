import sys
import time
import signal
import subprocess
import click
import pandas as pd

from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler, FileSystemEvent
from pathlib import Path

from live_sturgeon_plotting import plot_confidence_over_time, write_progress_tsv, plot_CNV_bam


# ==========================
# Constants
# ==========================
DEFAULT_WATCH_DIR: Path = Path("./bam_pass")
DEFAULT_SCRIPT: Path = Path("/Users/k.v.cammel/Developer/cold_setup_sturgeon/sturgeon_P2.sh")
DEFAULT_OUT_DIR: Path = Path("./sturgeon_results")
DEFAULT_SUFFIXES: str = ".bam"
DEFAULT_LOCK_FILE: Path = Path("script.lock")
DEFAULT_MODEL_FILE: Path = Path("/Users/k.v.cammel/Developer/cold_setup_sturgeon/sturgeon/include/models/general.zip")


# ==========================
# Lock File Utilities
# ==========================
def check_lock(lock_file: Path) -> None:
    """Check and create a lock file to prevent multiple instances."""
    if lock_file.exists():
        print("A different instance of Sturgeon is already running. Exiting...")
        sys.exit(1)
    with open(lock_file, "w") as lf:
        lf.write("Lock file to prevent multiple instances of sturgeon.")
    print("Lock file created. Proceeding with live processing...")


def remove_lock(lock_file: Path) -> None:
    """Remove the lock file on shutdown."""
    if lock_file.exists():
        lock_file.unlink()
        print("Lock file removed.")

def handle_exit(signum, frame):
    print("Received termination notice, exiting cleanly...")
    sys.exit(0)
    running = False
signal.signal(signal.SIGTERM, handle_exit)
signal.signal(signal.SIGINT, handle_exit)
# ==========================
# File Event Handler
# ==========================
class NewBamFileHandler(FileSystemEventHandler):
    def __init__(self, script_path: Path, results: Path, freq: int, model: Path) -> None:
        self.iteration = 1
        self.script_path = script_path
        self.output = results
        self.freq = freq
        self.model = model
        self.full_data = pd.DataFrame()

    def on_created(self, event: FileSystemEvent) -> None:
        """React to newly created files that match allowed suffixes."""
        if not event.is_directory:
            new_file = Path(event.src_path)
            filename = new_file.name
            if filename.endswith(DEFAULT_SUFFIXES):
                print(f"Detected new file: {filename}")
                self.run_script(new_file)

    def run_script(self, new_file: Path) -> None:
        """Run sturgeon and plotting on new file"""
        try:
            print(f"FLAG: starting processing of iteration_{self.iteration}")
            subprocess.run(["sh", self.script_path, new_file, self.output, self.model, str(self.iteration)], check=True)
            ##Copy output files with iteration number to out file
            if (self.iteration % self.freq == 0):
                self.plot_cnv(new_file)
            self.plot_process()
            print(f"FLAG: iteration_{self.iteration} completed!")
            self.iteration += 1
        except subprocess.CalledProcessError as e:
            print(f"Script failed with error: {e}")

    def plot_process(self) -> None:
        """Generate confidence over time plot and tsv"""
        print(f"FLAG: Making confidence over time plot for iteration_{self.iteration}", flush=True)
        modelname = self.model.stem
        self.full_data = write_progress_tsv(self.full_data, self.output, self.iteration, modelname)
        self.full_data.to_csv(self.output / f"iteration_{self.iteration}/classifier_progress_iteration_{self.iteration}.tsv",sep="\t",index=False)
        output_file = f"{self.output}/iteration_{self.iteration}/confidence_over_time_plot_iteration_{self.iteration}"
        df = pd.read_csv("/Users/k.v.cammel/Documents/Projects/P2_wrapper/color_translation.csv")
        color_translation = dict(zip(df["class"], df["color"]))
        plot_confidence_over_time(self.full_data, output_file, color_translation)



    def plot_cnv(self,new_file) -> None:
        """Generate CNV plot"""
        print(f"FLAG: Creating CNV plot for iteration {self.iteration}")
        ## R cannot handle Path type, so convert to string
        str_output = str(self.output)
        str_bam = new_file.absolute().as_posix()
        output_file = f"{str_output}/iteration_{self.iteration}/CNV_plot_iteration_{self.iteration}"
        plot_CNV_bam(str_bam,output_file)


# ==========================
# CLI with Click
# ==========================
@click.command()
@click.option(
    "-i",
    "--input",
    type=click.Path(path_type=Path, exists=False, dir_okay=True),
    default=DEFAULT_WATCH_DIR,
    help="Directory to watch for new sequencing result files.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path, exists=False, writable=True),
    default=DEFAULT_OUT_DIR,
    help="Directory where results are written.",
)
@click.option(
    "-l",
    "--lock",
    type=click.Path(path_type=Path, exists=False, writable=True),
    default=DEFAULT_LOCK_FILE,
    help="Name of lock file.",
)
@click.option(
    "-s",
    "--script",
    type=click.Path(path_type=Path, exists=True, file_okay=True),
    default=DEFAULT_SCRIPT,
    help="Path to the script that will be called for processing.",
)
@click.option(
    "-b",
    "--barcode",
    type=int,
    default=18,
    help="Barcode used in library preparation."
)
@click.option(
    "-f",
    "--freq",
    type=int,
    default=1,
    help="Number of iterations before merging BAMs and plotting CNV.",
)
@click.option(
    "-m",
    "--model",
    type=click.Path(path_type=Path, exists=True),
    default=DEFAULT_MODEL_FILE,
    help="Location of model used for sturgeon prediction"
)
def main(input: Path, output: Path, lock: Path, script: Path, barcode: int, freq: int, model: Path):
    """
    Monitors a folder for new BAM or TXT files and processes them as they appear.
    """
    check_lock(lock)
    while running:
        try:
            print(f"Starting to monitor for new BAM files in: {input}")
            event_handler = NewBamFileHandler(script_path=script, results=output, freq=freq, model=model)
            observer = Observer()
            observer.schedule(event_handler, path=input, recursive=False)
            observer.start()

            try:
                while True:
                    time.sleep(1)
            except KeyboardInterrupt:
                print("\nLive processing stopping...")
            finally:
                observer.stop()
                observer.join()
        finally:
            remove_lock(lock)


# ==========================
# Entry Point
# ==========================
if __name__ == "__main__":
    running = True
    main()
