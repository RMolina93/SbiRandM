import argparse, warnings
import .superimp as sup


def parse_arguments():

    """
    Parsing of arguments. 
    -d : Path to the folder with the PDB's Pairwise Interactions that conforms the Complex.
    -fasta: Path to the Fasta File you want of the Complex that you want to model.
    -output: Path to the folder where generate the output file complex.
    -v: Activate the verbosity mode
    """

    parser = argparse.ArgumentParser(description='SbiRandM Complex Builder 1.0.0. This program generates a Protein-DNA or Protein-Protein complex, based on a Fasta Sequence and the pairwise interaction structures of the complex.This version uses a Superimposition-based approach to generate the structure. More information can be found at https://crowthebio.gitbook.io/complex-builder')

    parser.add_argument('-d', action="store", required =True, dest="folder", help="Path to the folder with the PDB's Pairwise Interactions that conforms the Complex.")
    parser.add_argument('-fasta', action="store", required=True, dest="fasta_seq", help="Path to the Fasta File you want of the Complex that you want to model.")
    parser.add_argument('-output', action="store", required=True, dest="output_folder", help="Path to the folder where generate the output file complex.")
    parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False, help = 'Activate the verbosity mode')

    args = vars(parser.parse_args())
    return args

args = parse_arguments()

sup.mainSuperimp(args)
