import os, sys, warnings
import difflib, glob
from Bio.PDB import PDBParser
from .data import aminoacids
from collections import defaultdict
from .exceptions import *

def check_homology(fasta_1, fasta_2):
    """
    This function check the identity percentage between two chains, and returns True if they are homologous
    or False if they are not. The threshold is currently at 70%.

    @ Input - Fasta_1 : String of fasta sequence
              Fasta_2 : String of fasta sequence

    @ Output - Boolean variable indicating if sequence are homologous
    """

    if fasta_1 in fasta_2:
        return True
    elif fasta_2 in fasta_1:
        return True
    else:
        if difflib.SequenceMatcher(None, fasta_1, fasta_2).ratio() > 0.60:
            return True
        else:
            return False

def Get_fasta(chain):
    """
    This function takes as input a Biopython PDB Chain object, and parse it into 
    a sequence of aminoacids.

    @ Input - Biopython PDB Chain object
    @ Output - String of the sequence of Aminoacids of the PDB Chain Object.
    """

    sequence = ""
    for residue in chain:
        try:
            sequence += aminoacids[residue.get_resname().strip()]
        except:
            continue
    return sequence

def check_clash(structure_complex, mobile_chain):
    """
    This function checks if the new added chain clash with any chain of the complex.

    @ Input - Structure_complex : Biopython object of the actual complex
              Mobile_chain : Biopython object of the moving chain when adding a new chain.

    @ Output - True : Chain is clashing with complex.
               False: Chain is not clashing with complex. 
    """
    #return True 

    backbone_atoms = ["CB", "P"]
    for atom in list(structure_complex.get_atoms()):
        for atom_2 in list(mobile_chain.get_atoms()):
            if atom.id in backbone_atoms and atom_2.id in backbone_atoms:
                if atom_2 - atom < 1:
                    print ("Ups there is a crash!")
                    return False
    return True

def obtain_pairwise_dict(steichiometry_dict, TMP_folder):

    """
    This function parse a folder full with PDB of the pairwise interactions, and 
    parse them into a dictionary that contains as keys the interacting chains, and
    as value, the path of the file that contains the structure of the interaction.

    @Input - Steichiometry dict with the next format.
              Name of the chain : 
                        Steichiometry - Absolute number that the chain appear in the Fasta.
                        Sequence - Sequence of the chain

    @Output - Dictionary with the next keys:
              [Name of Chain] - [Name of chain in the interaction] - [Path to that PDB file]
    """

    pairwise_dict = dict()
    pairwise_dict = defaultdict(dict)

    parser = PDBParser(PERMISSIVE=1)
    index = 0

    for pdb_file in glob.glob(os.path.join(TMP_folder,"*.pdb")):
        index += 1
        structure = parser.get_structure('Complex', pdb_file)

        for steichiometry_chain in steichiometry_dict: #A,C
            for chain in structure.get_chains():
                if check_homology(Get_fasta(chain), steichiometry_dict[steichiometry_chain]["sequence"]):
                    chain.real_id = steichiometry_chain

        chain_list = list(structure.get_chains())

        try:
            pairwise_dict[chain_list[0].real_id][chain_list[1].real_id] = (pdb_file)
            pairwise_dict[chain_list[1].real_id][chain_list[0].real_id] = (pdb_file)

        except Exception:
            raise FilesDontMatchException

    return dict(pairwise_dict)

            
