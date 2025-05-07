import os
import time
from pathlib import Path


# === Config ===
SOURCE_DIR = Path("/Users/k.v.cammel/Documents/Projects/Sturgeon/P2_test_runs_250417/20250417_CNS_test_p3/no_sample_id/20250417_1421_P2I-00332-A_PBC28823_8dc9d25f/bam_pass/") ##Change this to own location, from prom@p22 /data/20250417_CNS_test_p3
DEST_DIR = Path("./bam_pass")
SYMLINK_PREFIX = "symlink_"
SYMLINK_INTERVAL = 2 * 60  # 10 minutes in seconds

last_index = -1

def create_symlink(source_file: Path, dest_folder: Path):
    """Create a symbolic link in the destination folder."""
    symlink_name = f"{SYMLINK_PREFIX}{source_file.name}"
    dest_file = dest_folder / symlink_name
    if dest_file.exists():
        print(f"Symlink {symlink_name} already exists. Skipping.")
    else:
        os.symlink(source_file, dest_file)
        print(f"Created symlink: {dest_file}")

def get_sorted_bam_files(source_dir: Path):
    """Get all BAM files in the source directory sorted alphabetically."""
    return sorted(source_dir.glob("*.bam"))

def main():
    # Ensure destination folder exists
    global last_index
    DEST_DIR.mkdir(parents=True, exist_ok=True)
    time.sleep(20)

    while True:
        # Get sorted BAM files from the source directory
        bam_files = get_sorted_bam_files(SOURCE_DIR)

        if bam_files:
            # Select the next BAM file (you can choose your strategy, e.g. the first, random, etc.)
            last_index = (last_index + 1) % len(bam_files)
            selected_bam_file = bam_files[last_index]  # You can randomize this if you want, e.g. random.choice(bam_files)
            print(f"Selected BAM file: {selected_bam_file}")

            # Create the symlink
            create_symlink(selected_bam_file, DEST_DIR)
        else:
            print("No BAM files found in source directory.")

        # Wait for the next cycle (2 minutes)
        print(f"Waiting for {SYMLINK_INTERVAL / 60} minutes...\n")
        time.sleep(SYMLINK_INTERVAL)

if __name__ == "__main__":
    main()
