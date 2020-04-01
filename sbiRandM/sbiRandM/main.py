import argparse, warnings
from sbiRandM.sbiRandM import *

warnings.filterwarnings("ignore")

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

#args = parse_arguments()

def mainSuperimp(args):
    if not os.path.isdir(args['output_folder']):
        os.mkdir(args['output_folder'])
    steichiometry_dict = check_fasta_stoichometry(args['fasta_seq'])
    pairwise_dict = obtain_pairwise_dict(steichiometry_dict , args['folder'])
    # print (pairwise_dict)
    pdb_complex = execute_complex(steichiometry_dict , pairwise_dict , args['output_folder'])
    
    while not model_validation(pdb_complex , args['fasta_seq']):
        print("Something went wrong, redoing model again.")
        try:
            pdb_complex = execute_complex(steichiometry_dict , pairwise_dict , args['output_folder'])
        except:
            print("Wow, something went VERY wrong! Redoing model")
            continue


#mainSuperimp(args)
