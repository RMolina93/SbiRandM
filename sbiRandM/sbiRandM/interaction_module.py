import os, sys, warnings
import difflib, glob
from Bio.PDB import PDBParser
from sbiRandM.sbiRandM import *
from sbiRandM.sbiRandM.data import aminoacids
from collections import defaultdict

def Get_fasta(chain):
    """
    @ Input - Biopython PDB Chain object
    @ Output - Sequence of Aminoacids of the PDB Chain Object.
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
    Check if the new added chain clash with the total complex.

    @ Input - Structure_complex : Biopython object of the actual complex
              Mobile_chain : Biopython object of the moving chain when adding a new chain.

    @ Output - True : Chain is clashing with complex.
               False: Chain is not clashing with complex. 
    """
    backbone_atoms = ["CB", "P"]
    for atom in list(structure_complex.get_atoms()):
        for atom_2 in list(mobile_chain.get_atoms()):
            if atom.id in backbone_atoms and atom_2.id in backbone_atoms:
                if atom_2 - atom < 1:
                    #print ("Ups there is a crash!")
                    return False
    return True

def obtain_pairwise_dict(steichiometry_dict, TMP_folder):

    """
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

        #print ("Parsing PDB file", pdb_file)

        for steichiometry_chain in steichiometry_dict: #A,C
            for chain in structure.get_chains():
                if difflib.SequenceMatcher(None, Get_fasta(chain), steichiometry_dict[steichiometry_chain]["sequence"]).ratio() > 0.70:
                    chain.real_id = steichiometry_chain
                    print("hola")
                print(steichiometry_chain)


        chain_list = list(structure.get_chains())

        pairwise_dict[chain_list[0].real_id][chain_list[1].real_id] = (pdb_file)
        pairwise_dict[chain_list[1].real_id][chain_list[0].real_id] = (pdb_file)

    print ("Parsed", index, "interaction files.")
    return dict(pairwise_dict)

            
