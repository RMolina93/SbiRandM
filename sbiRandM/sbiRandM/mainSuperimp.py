import argparse, warnings
from sbiRandM.sbiRandM import *
from sbiRandM.sbiRandM import superimp as sup


def parse_arguments():

    """
    Parsing of arguments. 
    -d : Folder with PDB's Pairwise Interactions
    -fasta: File with the fasta you want to model
    -output: Folder to put the output
    """

    parser = argparse.ArgumentParser(description='Complexes modelling application.')

    parser.add_argument('-d', action="store", required =True, dest="folder", help="Folder with all the PDBs of the interactions complex")
    parser.add_argument('-fasta', action="store", required=True, dest="fasta_seq", help="Fasta sequence of the full complex.")
    parser.add_argument('-output', action="store", required=True, dest="output_folder", help="Output dir to store the results")

    args = vars(parser.parse_args())
    return args

args = parse_arguments()

sup.mainSuperimp(args)
