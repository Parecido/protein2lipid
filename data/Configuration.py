import json
import os


class Configuration:
    def __init__(self, file_input, **kwargs):
        self.__dict__.update(kwargs)

        with open(f"build/{file_input}", "r") as file_in:
            input_json = json.load(file_in)

        for option in ["gromacs", "bilayer", "peptides"]:
            if option not in input_json:
                print(f"Missing {option} setup in input file. Please validate input file.")
                exit()

        for key, value in input_json.items():
            if key == "gromacs":
                self.mdrun = value["mdrun"]
                self.temperature = value["temperature"]
                self.pressure = value["pressure"]
            if key == "bilayer":
                self.template = value["template"]
                self.toppar = value["toppar"]
                self.lipids = [lipid for lipid in value["lipids"].split(",") if lipid.strip() != ""]
                self.additives = [additive for additive in value["additives"].split(",") if additive.strip() != ""]
            if key == "peptides":
                self.peptides = [peptide for peptide in value.values()]

        if len(self.lipids) == 0:
            print("Missing lipid names. Please validate input file.")
            exit()

        for lipid in self.lipids:
            if not os.path.exists(f"build/{self.toppar}/{lipid}.itp"):
                print(f"Missing {lipid}.itp file in {self.toppar} folder.")
                exit()

        for additive in self.additives:
            if not os.path.exists(f"build/{self.toppar}/{additive}.itp"):
                print(f"Missing {additive}.itp file in {self.toppar} folder.")
                exit()
