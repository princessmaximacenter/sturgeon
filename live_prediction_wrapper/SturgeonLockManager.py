import sys
from pathlib import Path


from live_prediction_wrapper import SturgeonLogging as SL


app_log = SL._get_app_logger()

class LockManager:
    def __init__(self, lock_file: Path):
        self.lock_file = lock_file

    def _check_lock(self) -> None:
        """
        Checks and creates lock file to prevent multiple instances
        :return: None
        """
        if self.lock_file.exists():
            app_log.error("A different instance of Sturgeon is already running. Exiting...")
            sys.exit(1)
        with open(self.lock_file, "w") as lf:
            lf.write("Lock file to prevent multiple instances of sturgeon.")
        app_log.info("Lock file created. Proceeding with live processing...")

    def _remove_lock(self) -> None:
        """
        Removes the lock file on shutdown
        :return: None
        """
        if self.lock_file.exists():
            self.lock_file.unlink()
            app_log.info("Lock file removed.")

