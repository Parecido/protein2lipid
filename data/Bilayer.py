from collections import OrderedDict
from pathlib import Path
import MDAnalysis as mda
import MDAnalysis.analysis.leaflet as lf
import os


class Bilayer:
    def __init__(self, filename, lipids, additives):
        self.file = open(f"build/{filename}", "r")
        self.lipids = lipids
        self.additives = additives
        self.universe = mda.Universe(self.get_absolute_path())

    def get_file(self):
        return Path(self.file.name)

    def get_base_name(self):
        return os.path.basename(self.get_file())

    def get_absolute_path(self):
        return self.get_file().resolve()

    def get_bilayer(self):
        return " ".join(self.lipids + self.additives)

    def get_atoms(self):
        return self.universe.select_atoms(f"resname {self.get_bilayer()}")

    def get_dimensions(self):
        return self.universe.dimensions

    def get_pattern(self):
        bilayer = self.lipids + self.additives
        residues = self.universe.residues
        return OrderedDict.fromkeys([residue.resname for residue in residues if residue.resname in bilayer])

    def get_residues(self):
        bilayer = dict()

        for lipid in self.lipids:
            atoms = self.universe.select_atoms(f"resname {lipid}")
            bilayer[lipid] = atoms

        for additive in self.additives:
            atoms = self.universe.select_atoms(f"resname {additive}")
            bilayer[additive] = atoms

        return list((residue, bilayer.get(residue)) for residue in list(self.get_pattern()))

    def get_surface(self, universe, to_convert):
        atom_p = "P31" if to_convert else "P"

        leaflet = lf.LeafletFinder(universe, f"resname {self.get_bilayer()} and name {atom_p}")
        leaflet0 = leaflet.groups(0)
        leaflet1 = leaflet.groups(1)

        z1 = leaflet0.atoms.center_of_mass()[2]
        z2 = leaflet1.atoms.center_of_mass()[2]

        return list(sorted([z1, z2]))
