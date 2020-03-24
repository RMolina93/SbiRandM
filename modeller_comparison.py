
# Homology modeling with multiple templates
from modeller import *              # Load standard Modeller classes
#from modeller.automodel import *    # Load the automodel class
import argparse
import glob, os, argparse, warnings
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from data import aminoacids

warnings.filterwarnings("ignore")
'''
parser = argparse.ArgumentParser(description='Complexes modelling application.')

parser.add_argument('-d', action="store", required =True, dest="folder", help="Folder with all the PDBs of the interactions complex")
parser.add_argument('-fasta', action="store", required=False, dest="fasta_seq", help="Fasta sequence of the full complex.")
parser.add_argument('-output', action="store", required=True, dest="output_folder", help="Output dir to store the results")

args = vars(parser.parse_args())
'''
def block_alignment(pir_file):
    """
    This function should return the alignment with the DNA BLOCKED for modeller.
    """
    return None


def obtain_sequence(folder):
    """
    This function takes a folder that contains the complexes of the pairwise interactions, 
    parse them and return a JSOn with information about the location and the sequences.

    @input folder - Folder where you store the pairwise complexes of the protein
    @output - JSON formatted variable with -  Proteins : chains : complex_path_file / sequence
    """

    templates_dict = dict()
    parser = PDBParser(PERMISSIVE=1)
    for pdb_file in glob.glob(os.path.join(folder,"*.pdb")):
        pdb_name = os.path.basename(pdb_file)
        templates_dict[pdb_name] = {"path":pdb_file, "name":pdb_name}
        structure = parser.get_structure('Complex', pdb_file)

        for pdb_chain in structure[0]:
            templates_dict[pdb_name][pdb_chain.id] = ""
            for residue in pdb_chain:
                if residue.get_resname() != 'UNK': # Por si en el pdb hay residuos desconocidos. SOLUCION A ARREGLAR
                    templates_dict[pdb_name][pdb_chain.id] += aminoacids[residue.get_resname()]
    print(templates_dict)
    return templates_dict


#if __name__ == "__main__":
#    obtain_sequence(args['folder'])


