import sys
import subprocess
import pandas as pd
from pathlib import Path
import os
from watchdog.events import FileSystemEventHandler, FileSystemEvent
from queue import Queue
import threading
import re
import pysam


import SturgeonLivePlotting as SLP
import SturgeonLogging as SL


log = SL._get_logger()


class LockManager:
    def __init__(self, lock_file: Path):
        self.lock_file = lock_file

    def _check_lock(self) -> None:
        """
        Checks and creates lock file to prevent multiple instances
        :return: None
        """
        if self.lock_file.exists():
            log.error("A different instance of Sturgeon is already running. Exiting...")
            sys.exit(1)
        with open(self.lock_file, "w") as lf:
            lf.write("Lock file to prevent multiple instances of sturgeon.")
        log.info("Lock file created. Proceeding with live processing...")

    def _remove_lock(self) -> None:
        """
        Removes the lock file on shutdown
        :return: None
        """
        if self.lock_file.exists():
            self.lock_file.unlink()
            log.info("Lock file removed.")


class NewBamFileHandler(FileSystemEventHandler):
    def __init__(self, sturgeon_script_path: Path, output: Path, model: Path, freq: int, utils: Path, r_script_path: Path, input: Path, gridion: bool) -> None:
        self.iteration = 1
        self.script_path = sturgeon_script_path
        self.output = output
        self.model = model
        self.freq = freq
        self.full_data = pd.DataFrame()
        self.utils = utils
        self.r_script_path = r_script_path
        self.watch_directory = input
        self.gridion = gridion

        """Create queue and processing thread for bam files (in chronological order)"""
        self.file_queue = Queue()
        self.processing_thread = threading.Thread(target=self._process_queue,daemon=True)
        self.processing_thread.start()

        """Process gridion run or already present bam files in P2 run"""
        if gridion:
            self._process_gridion_run()
        else:
            self._process_existing_results()

    def _process_existing_results(self):
        """
        Checks the output directory for bam files that were already present before Sturgeon was started.
        :return: None
        """
        log.info("Checking existing bam files")

        bamFiles = []
        for bamFile in self.watch_directory.glob("*.bam"):
            creationTime = os.path.getmtime(bamFile)

            bamFiles.append((bamFile,creationTime))

        bamFiles.sort(key=lambda x:x[1]) #sort by oldest bam file first

        for bamFile, _ in bamFiles:
            self.file_queue.put(bamFile)


    def _process_gridion_run(self) -> None:
        """
        Collect all the bam files from the gridion run and create symlinks for processing
        :return: None
        """
        log.info("Processing files from gridion run")

        bamFiles = list(self.watch_directory.glob("guppy_output_it*[0-9]*.bam"))
        mainBamFiles = []
        for bam in bamFiles:
            if "unclassified" in bam.name:
                continue
            iteration = re.search(r"it(?:_iteration)?_(\d+)|it(\d+)", bam.name)
            if iteration:
                it_str = iteration.group(1) or iteration.group(2)
                it_num = int(it_str)
                mainBamFiles.append((bam, it_num))

        #Check if unclassified barcodes were used in gridion analysis
        if any("unclassified" in bam.name for bam in bamFiles):
            mergedBams = []
            for mainBam, it_num in sorted(mainBamFiles,key=lambda x: x[1]):
                if "_iteration_" in mainBam.name:
                    unclassified_bam = self.watch_directory / f"guppy_output_it_iteration_{it_num}_unclassified.bam"
                else:
                    unclassified_bam = self.watch_directory / f"guppy_output_it{it_num}_unclassified.bam"


                mergedBam = [str(mainBam),str(unclassified_bam)]
                Path(f"{self.output}/merged_bams").mkdir(exist_ok=True)
                mergedPath = f"{self.output}/merged_bams/merged_it{it_num}.bam"
                pysam.merge("-O", "BAM", "-o", mergedPath, *mergedBam)
                mergedBams.append((mergedPath,it_num))

            for mergedPath, _ in sorted(mergedBams,key=lambda x: x[1]):
                self.file_queue.put(mergedPath)
        else:
            for bamPath, it_num in sorted(mainBamFiles, key=lambda x: x[1]):
                symlink_target = Path(f"{self.output}/merged_bams/merged_it{it_num}.bam")
                symlink_target.parent.mkdir(parents=True, exist_ok=True)
                symlink_target.symlink_to(Path(bamPath))
                self.file_queue.put(bamPath)

    def _process_queue(self):
        while True:
            filePath = self.file_queue.get()
            self.run_script(filePath)

    def on_created(self, event: FileSystemEvent) -> None:
        """
        Process newly created/detected bam file and adds to queue
        :param event: Newly detected bam file
        :return: None
        """
        if not event.is_directory and event.src_path.endswith(".bam"):
            log.info(f"Detected new file: {event.src_path}")
            self.file_queue.put(Path(event.src_path))


    def run_script(self, new_file: Path) -> None:
        """
        Function that runs the sturgeon bash script and all plotting functions (confidence + CNV) for file first in queue
        :param new_file: bam file processed from the queue
        :return: None
        """
        try:
            log.info(f"FLAG: Starting processing of iteration_{self.iteration}")
            subprocess.run(
                ["bash", str(self.script_path), str(new_file), str(self.output), str(self.model), str(self.iteration)],
                check=True
            )
            if self.gridion == False:
                ##Create symlink for merging of bams for CNV file
                symlink_target = Path(f"{self.output}/merged_bams/bam_for_CNV_it{self.iteration}.bam")
                symlink_target.parent.mkdir(parents=True, exist_ok=True)
                symlink_target.symlink_to(new_file)
            else:
                pass
            if self.iteration % self.freq == 0:
                Path(f"{self.output}/merged_bams").mkdir(exist_ok=True)
                self.plot_cnv()
            self.plot_process()
            log.info(f"FLAG: iteration_{self.iteration} completed!")
            self.iteration += 1
            log.info(f"FLAG: Waiting for new bam file")
        except subprocess.CalledProcessError as e:
            log.error(f"Script failed with error: {e}")

    def plot_process(self) -> None:
        """
        Plots the confidence over time plot for the current iteration
        :return: None
        """
        log.info(f"FLAG: Creating confidence over time plot for iteration_{self.iteration}")
        modelname = self.model.stem
        self.full_data = SLP.write_progress_tsv(self.full_data, self.output, self.iteration, modelname)
        self.full_data.to_csv(
            self.output / f"iteration_{self.iteration}/classifier_progress_iteration_{self.iteration}.tsv",
            sep="\t", index=False
        )
        output_file = f"{self.output}/iteration_{self.iteration}/confidence_over_time_plot_iteration_{self.iteration}"
        color_translation = self._load_color_translation()
        SLP.plot_confidence_over_time(self.full_data, output_file, color_translation)


    def _load_color_translation(self) -> dict:
        df = pd.read_csv(f"{self.utils}/color_translation.csv") #Set in config.yaml
        return dict(zip(df["class"], df["color"]))

    def plot_cnv(self) -> None:
        """
        Merges bam files to create bam file used for plotting the CNV for the current iteration
        :return: None
        """
        log.info(f"FLAG: Creating CNV plot for iteration_{self.iteration}")
        bamToCNV = f"{self.output}/merged_bams/merged_CNV_bam.bam"
        bam_dir = Path(self.output) / "merged_bams"
        if self.gridion == True:
            modkitBamsToMerge = sorted(
                [
                    bam for bam in bam_dir.glob("merged_it*.bam")
                    if (match := re.search(r'merged_it(\d+)\.bam$', bam.name)) and int(match.group(1)) <= self.iteration
                ],
                key=lambda bam: int(re.search(r'merged_it(\d+)\.bam$', bam.name).group(1))
            )
        else:
            modkitBamsToMerge = sorted(
                [
                    bam for bam in bam_dir.glob("bam_for_CNV_it*.bam")
                    if (match := re.search(r'bam_for_CNV_it(\d+)\.bam$', bam.name)) and int(match.group(1)) <= self.iteration
                ],
                key=lambda bam: int(re.search(r'bam_for_CNV_it(\d+)\.bam$', bam.name).group(1))
            )

        pysam.merge("-@","10","-f","-O","BAM","-o", bamToCNV, *map(str, modkitBamsToMerge)) #Merge all bams after modkit, before creating CNV plot
        output_file = f"{self.output}/iteration_{self.iteration}/CNV_plot_iteration_{self.iteration}"
        SLP.plot_CNV_bam(bamToCNV, output_file, self.r_script_path,self.utils)




