import sys
import subprocess
import pandas as pd
from pathlib import Path
import os
from watchdog.events import FileSystemEventHandler, FileSystemEvent
from queue import Queue
import threading

import SturgeonLivePlotting as SLP
import SturgeonLogging as SL


log = SL._get_logger()


class LockManager:
    def __init__(self, lock_file: Path):
        self.lock_file = lock_file

    def check_lock(self) -> None:
        """Check and create a lock file to prevent multiple instances."""
        if self.lock_file.exists():
            log.error("A different instance of Sturgeon is already running. Exiting...")
            sys.exit(1)
        with open(self.lock_file, "w") as lf:
            lf.write("Lock file to prevent multiple instances of sturgeon.")
        log.info("Lock file created. Proceeding with live processing...")

    def remove_lock(self) -> None:
        """Remove the lock file on shutdown."""
        if self.lock_file.exists():
            self.lock_file.unlink()
            log.info("Lock file removed.")


class NewBamFileHandler(FileSystemEventHandler):
    def __init__(self, sturgeon_script_path: Path, output: Path, model: Path, freq: int, r_script_path: Path, input: Path) -> None:
        self.iteration = 1
        self.script_path = sturgeon_script_path
        self.output = output
        self.model = model
        self.freq = freq
        self.full_data = pd.DataFrame()
        self.r_script_path = r_script_path
        self.watch_directory = input

        """Create queue and processing thread for bam files (in chronological order)"""
        self.file_queue = Queue()
        self.processing_thread = threading.Thread(target=self._process_queue,daemon=True)
        self.processing_thread.start()

        """Process existing files"""
        self._process_existing_results()

    def _process_existing_results(self):
        """Check the results directory if there are already bam files present"""
        log.info("Checking existing bam files")

        bamFiles = []
        for bamFile in self.watch_directory.rglob("*.bam"):
            creationTime = os.path.getmtime(bamFile)

            bamFiles.append((bamFile,creationTime))

        bamFiles.sort(key=lambda x:x[1]) #sort by oldest bam file first

        for bamFile, _ in bamFiles:
            self.file_queue.put(bamFile)

    def _process_queue(self):
        while True:
            filePath = self.file_queue.get()
            self.run_script(filePath)

    def on_created(self, event: FileSystemEvent):
        """React to newly created files that match allowed suffixes."""
        if not event.is_directory and event.src_path.endswith(".bam"):
            log.info(f"Detected new file: {event.src_path}")
            self.file_queue.put(Path(event.src_path))


    def run_script(self, new_file: Path) -> None:
        """Run sturgeon and plotting on new file"""
        try:
            log.info(f"FLAG: Starting processing of iteration_{self.iteration}")
            subprocess.run(
                ["bash", str(self.script_path), str(new_file), str(self.output), str(self.model), str(self.iteration)],
                check=True
            )
            if self.iteration % self.freq == 0:
                self.plot_cnv(new_file)
            self.plot_process()
            log.info(f"FLAG: iteration_{self.iteration} completed!")
            self.iteration += 1
            log.info(f"FLAG: Waiting for new bam file")
        except subprocess.CalledProcessError as e:
            log.error(f"Script failed with error: {e}")

    def plot_process(self) -> None:
        """Generate confidence over time plot and tsv"""
        log.info(f"FLAG: Creating confidence over time plot for iteration_{self.iteration}")
        modelname = self.model.stem
        self.full_data = SLP.write_progress_tsv(self.full_data, self.output, self.iteration, modelname)
        self.full_data.to_csv(
            self.output / f"iteration_{self.iteration}/classifier_progress_iteration_{self.iteration}.tsv",
            sep="\t", index=False
        )
        output_file = f"{self.output}/iteration_{self.iteration}/confidence_over_time_plot_iteration_{self.iteration}"
        color_translation = self.load_color_translation()
        SLP.plot_confidence_over_time(self.full_data, output_file, color_translation)


    def load_color_translation(self):
        df = pd.read_csv("/Users/k.v.cammel/Developer/cold_setup_sturgeon/utils/color_translation.csv") #Set in config.yaml
        return dict(zip(df["class"], df["color"]))

    def plot_cnv(self, new_file) -> None:
        """Generate CNV plot"""
        log.info(f"FLAG: Creating CNV plot for iteration_{self.iteration}")
        str_output = str(self.output)
        str_bam = new_file.absolute().as_posix()
        output_file = f"{str_output}/iteration_{self.iteration}/CNV_plot_iteration_{self.iteration}"
        SLP.plot_CNV_bam(str_bam, output_file, self.r_script_path)




