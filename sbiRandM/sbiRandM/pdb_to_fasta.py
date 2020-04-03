import os, sys
from Bio.PDB import PDBParser
from sbiRandM.sbiRandM import *
from sbiRandM.sbiRandM.data import aminoacids, Alphabet

def pdb_to_fasta(pdb_file, out_fasta_file):
    """
    This function takes as input a PDB file, and
    returns the Fasta file of the sequence.
    """

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('Complex', pdb_file)

    index = 0
    with open(out_fasta_file, "w") as out_fasta:
        for chain in structure.get_chains():
            out_fasta.write(">" + str(Alphabet[index]) + "\n")
            out_fasta.write(Get_fasta(chain) + "\n")
            index += 1
    
    return out_fasta_file


