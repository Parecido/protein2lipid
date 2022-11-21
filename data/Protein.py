from pathlib import Path
import os


class Protein:
    def __init__(self, filename):
        self.file = open(f"build/{filename}", "r")

    def get_file(self):
        return Path(self.file.name)

    def get_base_name(self):
        return os.path.basename(self.get_file())

    def get_stem(self):
        return Path(self.get_file()).stem
