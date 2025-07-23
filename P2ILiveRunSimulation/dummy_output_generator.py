import os
import time
from pathlib import Path
import yaml
import click

# Load config file with default values
pythonPath = Path(__file__).resolve()
TestingRoot = pythonPath.parents[0]
CONFIG_PATH = TestingRoot / "Config.yaml"

def load_config(config_yaml: Path) -> dict:
    with open(config_yaml, 'r') as config_file:
        return yaml.safe_load(config_file)

def create_symlink(source_file: Path, dest_folder: Path, prefix:str):
    """Create a symbolic link in the destination folder."""
    symlink_name = f"{prefix}{source_file.name}"
    dest_file = dest_folder / symlink_name
    if dest_file.exists():
        print(f"Symlink {symlink_name} already exists. Skipping.")
    else:
        os.symlink(source_file, dest_file)
        print(f"Created symlink: {dest_file}")

def get_sorted_bam_files(source_dir: Path):
    """Get all BAM files in the source directory sorted alphabetically."""
    return sorted(source_dir.glob("*.bam"))

@click.command()
@click.option('--source-dir', type=click.Path(exists=True, file_okay=False, path_type=Path), help="Source directory containing BAM files.")
@click.option('--dest-dir', type=click.Path(file_okay=False, path_type=Path), help="Destination directory where symlinks will be created.")
@click.option('--interval', type=int, default=60, show_default=True, help="Interval in seconds between each symlink creation.")
@click.option('--prefix', type=str, default="symlink_", show_default=True, help="Prefix given to symlinks")
@click.option('--config', type=click.Path(exists=True, path_type=Path), default=CONFIG_PATH, show_default=True, help="Path to YAML config file.")

def main(source_dir: Path, dest_dir: Path, interval: int, prefix:str, config: Path):
    # Load config from YAML
    config_values = load_config(config)

    # Resolve values: CLI > config > fail if required
    source_dir = source_dir or Path(config_values.get("source_dir", ""))
    dest_dir = dest_dir or Path(config_values.get("dest_dir", ""))
    interval = interval or config_values.get("interval", 60)
    prefix = prefix or config_values.get("prefix", "symlink_")
    # Ensure destination folder exists
    last_index = -1
    dest_dir.mkdir(parents=True, exist_ok=True)
    time.sleep(10)


    while True:
        # Get sorted BAM files from the source directory
        bam_files = get_sorted_bam_files(source_dir)

        if bam_files:
            # Select the next BAM file
            last_index = (last_index + 1) % len(bam_files)
            selected_bam_file = bam_files[last_index]
            print(f"Selected BAM file: {selected_bam_file}")

            # Create the symlink
            create_symlink(selected_bam_file, dest_dir,prefix)
        else:
            print("No BAM files found in source directory.")

        # Wait for the next cycle (2 minutes)
        print(f"Waiting for {interval / 60} minutes...\n")
        time.sleep(interval)

if __name__ == "__main__":
    main()
