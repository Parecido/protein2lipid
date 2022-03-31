from data.lipid2lipid import convert_bilayer
from MDAnalysis.analysis import align
import gromacs as gmx
import MDAnalysis as mda
import os
import shutil


def check_lines(line, skip_conditions):
    for skip_condition in skip_conditions:
        if skip_condition in line:
            return True, skip_condition

    return False, ""


def rewrite_file(file_in, file_out, skip_conditions):
    write_lines = []
    with open(file_in, "r") as file_r:
        for line in file_r:
            to_skip, key = check_lines(line, skip_conditions)
            if to_skip:
                next_line = next(file_r, "")
                while skip_conditions[key] != next_line.strip():
                    next_line = next(file_r, "")
            else:
                write_lines.append(line)

    with open(file_out, "w") as file_w:
        file_w.writelines(write_lines)


def prepare_dirs(protein, bilayer, conf):
    if os.path.isdir(f"build/{protein.get_stem()}"):
        shutil.rmtree(f"build/{protein.get_stem()}")

    if os.path.isdir("build/scratch"):
        shutil.rmtree("build/scratch")

    os.system("cp -r tools/ions build/scratch")
    os.system(f"cp {protein.get_file()} build/scratch/")
    os.system(f"cp {bilayer.get_file()} build/scratch/")
    os.system("cp -r tools/charmm36-jul2021.ff build/scratch/")
    os.system(f"cp -r build/{conf.toppar} build/scratch/")
    os.chdir("build/scratch")


def prepare_protein_topology(protein, conf):
    gmx.editconf(f=protein.get_base_name(), o=protein.get_base_name(), resnr=1)
    gmx.pdb2gmx(f=protein.get_base_name(), o="1protein.gro", ignh=True, merge="all",
                renum=True, inter=conf.prot_stat, ter=conf.terminal, input=conf.stdin)

    rewrite_file("topol.top", "1protein.top", {"[ system ]": "", "[ molecules ]": ""})
    return "1protein.gro", "1protein.top"


def aligment_protein(protein_gro, bilayer):
    mobile = mda.Universe(protein_gro)
    ref = mda.Universe(bilayer.get_base_name())

    align.alignto(mobile, ref, select=f"resid 1-{len(mobile.residues)} and name CA")
    return mobile.atoms


def merge_gros(protein_atoms, bilayer, graphs):
    bilayer_atoms = bilayer.get_atoms()

    if len(graphs) > 0:
        bilayer_atoms = convert_bilayer(bilayer, graphs)

    protein_bilayer = mda.Merge(protein_atoms, bilayer_atoms)
    protein_bilayer.dimensions = bilayer.get_dimensions()
    protein_bilayer.atoms.write("2merged.gro")

    gmx.editconf(f="2merged.gro", o="2merged-fixed.gro", resnr=1)
    return "2merged-fixed.gro"


def prepare_topology(protein_top, bilayer, conf):
    with open(protein_top, "r") as file_top:
        top_gromacs = file_top.readlines()

    lipids = bilayer.get_residues()

    if conf.to_convert:
        top_gromacs.insert(21, """#include "toppar/ffbonded.itp"\n""")
        os.system("cp ../../tools/converter/ffbonded.itp toppar/.")

        for lipid in lipids:
            os.system(f"cp ../../tools/converter/{lipid[0]}.itp toppar/.")

    top_gromacs += [f"""#include "toppar/{lipid[0]}.itp"\n""" for lipid in lipids]
    top_gromacs.append("\n[ system ]\n")
    top_gromacs.append("; Name\n")
    top_gromacs.append("Title\n")
    top_gromacs.append("\n")
    top_gromacs.append("[ molecules ]\n")
    top_gromacs.append("; Compound           #mols\n")
    top_gromacs.append("Protein_chain_A 1\n")
    top_gromacs += [f"{lipid[0]} {len(lipid[1].residues)}\n" for lipid in lipids]

    with open("2merged-fixed.top", "w") as f_merged:
        f_merged.writelines(top_gromacs)

    return "2merged-fixed.top"


def add_solvent(merged_gro, merged_top, bilayer, conf):
    gmx.solvate(cp=merged_gro, cs="spc216.gro", p=merged_top, o="3system_water.gro")

    universe = mda.Universe("3system_water.gro")
    z1, z2 = bilayer.get_surface(universe, conf.to_convert)
    water_bilayer = universe.select_atoms(f"(prop z >= {z1} and prop z <= {z2}) and resname SOL")

    water_removed = universe.select_atoms("not byres (group water_bilayer)", water_bilayer=water_bilayer)
    water_removed.atoms.write("3system_water.gro")

    with open(merged_top, "r") as f_top:
        top_water = f_top.readlines()

    water_atoms = water_removed.select_atoms("resname SOL")
    water_molecules = int(len(water_atoms) / 3)
    top_water = top_water[0:-1] + [f"SOL {water_molecules}\n"]

    with open("3system_water.top", "w") as f_water:
        f_water.writelines(top_water)

    return "3system_water.gro", "3system_water.top"


def add_ions(water_gro, water_top, conf):
    gmx.grompp(f="generic.mdp", c=water_gro, p=water_top, o="4ions")

    if conf.to_convert:
        gmx.genion(s="4ions.tpr", p=water_top, o="4system_ions.gro", conc=0.150, neutral=True, input="SOL")
    else:
        gmx.genion(s="4ions.tpr", p=water_top, o="4system_ions.gro", conc=0.150, neutral=True, input="SOL",
                   pname="SOD", nname="CLA")

    return "4system_ions.gro", "3system_water.top"


def prepare_output(system_gro, system_top, protein, bilayer, conf):
    os.system(f"cp -r ../../tools/mdp ../{protein.get_stem()}")
    os.system(f"cp {system_gro} ../{protein.get_stem()}/step5_input.gro")
    os.system(f"cp {system_top} ../{protein.get_stem()}/topol.top")
    os.system(f"cp posre.itp ../{protein.get_stem()}/")
    os.system(f"mkdir ../{protein.get_stem()}/toppar")

    if not conf.to_convert:
        os.system(f"cp -r charmm36-jul2021.ff ../{protein.get_stem()}/")
    else:
        os.system(f"cp toppar/ffbonded.itp ../{protein.get_stem()}/toppar/.")

    lipids = bilayer.get_residues()
    for lipid in lipids:
        file_itp = f"../{protein.get_stem()}/toppar/{lipid[0]}.itp"
        os.system(f"cp toppar/{lipid[0]}.itp ../{protein.get_stem()}/toppar/")
        rewrite_file(file_itp, f"{file_itp}_new", {"#ifdef": "#endif"})
        os.system(f"rm {file_itp}")
        os.system(f"mv {file_itp}_new {file_itp}")

    for file in os.listdir(f"../{protein.get_stem()}/"):
        if file == "README":
            file_readme = f"../{protein.get_stem()}/README"
            if conf.mdrun.strip() != "":
                os.system(f"sed -i 's/mdrun/mdrun {conf.mdrun}/g' {file_readme}")

        if file.endswith(".mdp"):
            file_mdp = f"../{protein.get_stem()}/{file}"
            os.system(f"sed -i 's/999 999/{conf.temperature} {conf.temperature}/g' {file_mdp}")
            os.system(f"sed -i 's/9.9 9.9/{conf.pressure} {conf.pressure}/g' {file_mdp}")

    os.chdir("../")
    shutil.rmtree("scratch")
    os.chdir("../")
