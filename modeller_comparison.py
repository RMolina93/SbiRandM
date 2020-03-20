
from modeller import *              
from modeller.automodel import *    
import argparse
import glob, os, argparse, warnings
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from data import aminoacids
from Bio import SeqIO
from models import Protein_Interaction, Chain, Query
import difflib

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

def fasta_to_object(fasta):

    """
    This function takes a fasta file with several chains, and parse it into a
    query object. Also changes the DNA chains to Placeholder characters.

    @input folder - File with Fasta
    @output - Query object.
    """

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    query = Query(name=os.path.basename(fasta))
    for sequence in record_dict.keys():
        chain = Chain(name = sequence, sequence = record_dict[sequence].seq)
        #print (sequence, ":",record_dict[sequence].seq)
        chain.dna_to_placeholder()
        query.add_chain(chain)

        print (chain.name, chain.sequence)
    return query

def create_models(folder):

    """
    This function takes a folder that contains the complexes of the pairwise interactions, 
    parse them and return a list of Protein_Interaction objects with information 
    about the location and the sequences.

    @input folder - Folder where you store the pairwise complexes of the protein
    @output - List of Protein_Interaction objects
    """

    list_of_interactions = list()
    parser = PDBParser(PERMISSIVE=1)
    for pdb_file in glob.glob(os.path.join(folder,"*.pdb")):
        pdb_name = os.path.basename(pdb_file)
        protein = Protein_Interaction(name = pdb_name, path = pdb_file)
        structure = parser.get_structure('Complex', pdb_file)

        for pdb_chain in structure[0]:
            sequence = ""
            for residue in pdb_chain:
                sequence += aminoacids[residue.get_resname().strip()]
            protein.add_chain(Chain(pdb_chain.get_id(), sequence))

        list_of_interactions.append(protein)
    
    return list_of_interactions

def check_similarity(query, interactions_list):

    """
    This function takes a query object and a list of Pairwise interaction objects, and
    returns the pairs of chains that have a similarity ratio on sequence over 73%
    """
    for chain in query.chains:
        for protein in interactions_list:
            for protein_chain in protein.chains:
                #print (chain.sequence.strip(), "/////", str(protein_chain.sequence)) 
                if difflib.SequenceMatcher(None,chain.sequence.strip(), str(protein_chain.sequence)).ratio() >0.73 :
                    print ("These two chains match", chain.name, protein.name, protein_chain.name)

if __name__ == "__main__":

    args = parse_arguments()
    query = fasta_to_object(args['fasta_seq'])
    interactions = create_models(args['folder'])
    check_similarity(query, interactions)


