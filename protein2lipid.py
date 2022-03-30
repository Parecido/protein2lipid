from data.lipid2lipid import create_graphs
from data.Tools import prepare_dirs, prepare_protein_topology, aligment_protein, merge_gros
from data.Tools import prepare_topology, add_solvent, add_ions, prepare_output
from optparse import OptionParser
import data.Bilayer as bl
import data.Input as ip
import data.Protein as pt


def parse_proteins(file_input="input.json", prot_stat=False, terminal=False,
                   to_convert=False, skip_symmetry=False, stdin=""):

    conf = ip.Input(file_input)
    bilayer = bl.Bilayer(conf.template, conf.lipids, conf.additives)
    graphs = create_graphs(bilayer, skip_symmetry) if to_convert else dict()

    for peptide in conf.peptides:
        protein = pt.Protein(peptide)

        prepare_dirs(protein, bilayer, conf)
        protein_gro, protein_top = prepare_protein_topology(protein, prot_stat, terminal, stdin)
        protein_atoms = aligment_protein(protein_gro, bilayer)
        merged_gro = merge_gros(protein_atoms, bilayer, graphs)
        merged_top = prepare_topology(protein_top, bilayer, to_convert)
        water_gro, water_top = add_solvent(merged_gro, merged_top, bilayer, to_convert)
        ions_gro, ions_top = add_ions(water_gro, water_top, to_convert)
        prepare_output(ions_gro, ions_top, protein, bilayer, conf, to_convert)


def main():
    parser = OptionParser()
    parser.add_option("-f", "--file_input", dest="file_input",
                      help="Input file name (default input.json)")
    parser.add_option("-p", "--prot_stat", action="store_true", dest="prot_stat",
                      help="Assign protonation states (default false)")
    parser.add_option("-t", "--terminal", action="store_true", dest="terminal",
                      help="Perform termini modifications (default false)")
    parser.add_option("-c", "--convert", action="store_true", dest="to_convert",
                      help="Perform membrane topology conversion (default false)")
    parser.add_option("-s", "--skip_symmetry", action="store_true", dest="skip_symmetry",
                      help="Skip molecular symmetry in topology conversion (default false")
    parser.add_option("-o", "--stdin", dest="stdin",
                      help="Autoselect topology parameters")
    (options, args) = parser.parse_args()

    parse_proteins(file_input=options.file_input, prot_stat=options.prot_stat,
                   terminal=options.terminal, to_convert=options.to_convert,
                   skip_symmetry=options.skip_symmetry, stdin=options.stdin)


if __name__ == "__main__":
    main()
