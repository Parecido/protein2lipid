import MDAnalysis as mda
import networkx as nx
import numpy as np
import os


def create_graphs(bilayer, skip_symmetry):
    graphs = dict()

    for residue in bilayer.get_residues():
        residue_name, _ = residue

        if not os.path.exists(f"tools/converter/{residue_name}.itp"):
            print(f"Missing {residue_name}.itp file for {residue_name} residue.")
            exit()

        itp1 = mda.Universe(f"tools/converter/{residue_name}.itp")
        itp2 = mda.Universe(f"build/toppar/{residue_name}.itp")

        bonds1 = [[bond.atoms.ids[0], bond.atoms.ids[1]] for bond in itp1.bonds]
        bonds2 = [[bond.atoms.ids[0], bond.atoms.ids[1]] for bond in itp2.bonds]

        G1 = nx.Graph(bonds1)
        G2 = nx.Graph(bonds2)

        ismags = nx.isomorphism.ISMAGS(G1, G2)
        largest_common_subgraph = list(ismags.largest_common_subgraph())

        if len(largest_common_subgraph) > 1 and not skip_symmetry:
            print("Detected symmetry in molecular graph. Aborting force field conversion.")
            exit()

        atoms = [(key, value, itp1.atoms[key-1].name) for key, value in largest_common_subgraph[0].items()]
        graphs[residue_name] = sorted(atoms, key=lambda v: v[1])

    return graphs


def convert_atoms(converter, residues, universe):
    residue_converted = mda.AtomGroup([], universe)

    for residue in residues:
        atoms_origin = universe.select_atoms(f"resid {residue.resid}")
        atoms_converted = np.zeros([len(atoms_origin)], dtype=object)

        for i, atom in enumerate(atoms_origin):
            new_index, _, new_name = converter[i]
            atom.name = new_name
            atoms_converted[new_index-1] = atom

        residue_converted += mda.AtomGroup([atom for atom in atoms_converted])

    return residue_converted


def convert_bilayer(bilayer, graphs):
    atoms_group = mda.AtomGroup([], bilayer.universe)

    for residue in bilayer.get_residues():
        residue_name, residue_atoms = residue
        atoms_group += convert_atoms(graphs[residue_name], residue_atoms.residues, bilayer.universe)

    return atoms_group
