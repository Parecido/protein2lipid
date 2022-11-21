# protein2lipid

Python framework for membrane topology conversion and auto-generation of membrane-bounded protein systems.

## Getting Started

### Dependencies

* Python (Anaconda)
* OS: Linux, MacOS
* Libraries: Gromacs (2018+), GromacsWrapper, MDAnalysis, NetworkX, NumPy

### Installing

* Download, compile, and install Gromacs 2018+
* Download latest protein2lipid version
* Prepare conda environment
```
conda create -n protein2lipid --file=requirements.txt -c bioconda
```

### Required files

All files should be located within build directory.

#### Template file

To run protein2lipid you need template file which can be easly created by [CHARMM-GUI](https://charmm-gui.org/). Template file contains atomic positions and describes protein-membrane system. By default [CHARMM-GUI](https://charmm-gui.org/) name it as step5_input.gro.

#### Topology folder

Topology folder should contain .itp files for all lipids localised inside membrane. Our software supports toppar folder generated by [CHARMM-GUI](https://charmm-gui.org/). The name of topology folder must be specified in input file.

#### Protein structure files

To perform auto-generation you should provide protein structures, preferably in pdb file format. The protein2lipid will prepare system consisted of protein, membrane, water, ions (NaCl 150 mmol/L), and their topology. Peptides must have fewer or the same amino acids number as the peptide used in the template file.

#### Input file

Information about peptides, composition of lipid bilayer, temperature, pressure, etc. should be included in input file. An example of input file is shown below.
```
{
    "gromacs": {
        "mdrun": "",
        "temperature": "310",
        "pressure": "1.0"
    },
    "bilayer": {
        "template": "step5_input.gro",
        "toppar": "toppar",
        "lipids": "POPC,POPE",
        "additives": "CHL1"
    },
    "peptides": {
        "peptide1": "AGYLLPKINLKPLAKLPKKIL.pdb",
        "peptide2": "FFHHIFRGIVHVGKTIHRLVTG.pdb",
        "peptide3": "FFSLIPSLVGGLISAFK.pdb"
    }
}
```

### Executing program

#### Activate an anaconda environment

```
conda activate protein2lipid
```

#### Executing protein2lipid from the command line

```
python protein2lipid.py -f input.json
```

#### Executing protein2lipid from another Python script

```
from protein2lipid import parse_proteins

parse_proteins(file_input="input.json")
```

#### Features

The protein2lipid supports basic protein adjustments such as termini modifications and protonation states determination. The auto-generation of protein topology might be controlled by --stdin option.

In order to perform topology conversion (by default between CHARMM36 and Lipid14 force fields) use --convert. We included CHL1, DLPC, DMPC, DOPC, DPPC, DPPE, POPC, and POPE Lipid14 topology files. If you want to convert between other force fields, you should provide .itp files with the same name as used in template file and place it inside tools/converter folder.

## Help

For more information please see:
```
python protein2lipid.py -h
```

## Authors

Michal Michalski ([Git](https://github.com/Parecido), [LinkedIn](https://www.linkedin.com/in/michal-michalski95), [Orcid](https://orcid.org/0000-0001-6969-2074), [ResearchGate](https://www.researchgate.net/profile/Michal-Michalski-5))

## Version History

* Version 1.1
    * Update CHARMM36M force field files
    * Load all pdb files by default

* Version 1.0
    * Initial release
    * Support for auto-generation of multiple protein-membrane systems 
    * Support for N/C-terminal modifications
    * Support for determination of protonation states
    * Implement CHARMM36 -> Lipid14 (Amber) topology conversion

