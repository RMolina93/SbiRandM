import sys, os, warnings
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from Bio import SeqIO
from .exceptions import *
from .interaction_module import *
from .data import Alphabet
import itertools, random
from shutil import copyfile
import difflib
import warnings

warnings.filterwarnings("ignore")


def Compute_equal_chain(structure_1, structure_2):
    
    """
    @ Input - Two Biopython PDB structures of an interaction.
    @ Output - A list of Atoms of same size that correspond to the same chain.
    """

    for chain in structure_1[0]:
        for chain_2 in structure_2[0]:

            fasta_1 = Get_fasta(chain)
            fasta_2 = Get_fasta(chain_2)
            print ("Checking:\n",fasta_1,"\n and \n", fasta_2)
            if check_homology(Get_fasta(chain), Get_fasta(chain_2)):
                print ("They are homologous")

                if len(fasta_1) > len(fasta_2):
                    smaller = chain_2
                    larger = chain
                else:
                    smaller = chain
                    larger = chain_2

                s = difflib.SequenceMatcher(None, fasta_1, fasta_2)
                miau = s.get_matching_blocks()
                i = miau[0].a
                size = miau[0].size
                print ("init:", i, "size:",size)

                residues = list(larger.get_residues())[i : i + size]
                print (len(residues), len(Get_fasta(smaller)))

                atoms = list()
                for residue in residues:
                    atoms = atoms + list(residue.get_atoms())

                if smaller == chain_2:
                    return atoms, list(smaller.get_atoms())
                else:
                    return list(smaller.get_atoms()), atoms

    print ("They are not homologous :(")
    return (None, None)

def check_chain_addition(complex_pdb, chain, pairwise_dict, steichiometry_dict):
    """
    Check if complex PDB and Adding Chain have atoms in common
    Output - True / False
    """
    warnings.filterwarnings("ignore")

    parser = PDBParser(PERMISSIVE=1)

    complex_pdb = parser.get_structure('Complex', complex_pdb)
    
    try:
        chains_to_test_in_complex = list(pairwise_dict[chain].keys())
    except KeyError:
        raise FastaRaroException
    common_chains = list()
    
    for chain_pdb in complex_pdb.get_chains():
        fasta = Get_fasta(chain_pdb)
        for chain_2 in chains_to_test_in_complex:
            if check_homology(fasta, steichiometry_dict[chain_2]["sequence"]):
                common_chains.append(pairwise_dict[chain][chain_2])
    
    return common_chains

    return None

def build_complex(file_1, file_2):
    """
    This function takes the complex output file (or in the first iteration one of the pairwise interactions)
    and another pairwise interaction PDB complex. Then it tries to add the chain to the complex until there is not clash

    @ Input - Two file path for a PDB interactions.
    @ Output - File path of the complex PDB file / Error: Chain cannot be added.
    """

    parser = PDBParser(PERMISSIVE=1)

    structure_1 = parser.get_structure('Complex', file_1)
    structure_2 = parser.get_structure('Complex', file_2)

    sup = Superimposer()
    io = PDBIO()


    atoms_fixed, atoms_moving = Compute_equal_chain(structure_1, structure_2)

    try:
        sup.set_atoms(atoms_fixed, atoms_moving)
    except:
        return False

    sup.apply(list(structure_2.get_atoms()))

    for chain in structure_2[0].get_chains():
        if chain.id != list(atoms_moving)[0].get_full_id()[2]:
            moved_chain = chain

    if check_clash(structure_1, moved_chain):
        with open(file_1, "wt") as out_file:
            
            for model in list(structure_1.get_chains()) + [moved_chain] :
                io.set_structure(model)
                io.save(out_file)

        rename_complex_chains(file_1)

        return True
    return False


def rename_complex_chains(file):
    """
    This function takes a PDB file, delete all the END lines, and rename the chains.
    according to Alphabet found in data module.

    @ Input - PDB file
    @ Output - PDB file with chains named.
    """

    index = 0
    temp_file = os.path.join(os.path.dirname(file),"temporal.pdb")
    with open(file, "r") as input_file:
        with open(temp_file,"w") as out_file:

            for line in input_file:
                if "TER" in line:
                    index += 1
                    continue
                else:
                    if "ATOM" in line:
                        line = list(line)
                        line[21] = Alphabet[index]
                        out_file.write("".join(line))
    os.remove(file)
    os.rename(temp_file, file)

def execute_complex(steichiometry_dict, pairwise_dict, args):
    warnings.filterwarnings("ignore")

    chains_list = [[i] * steichiometry_dict[i]["steichiometry"] for i  in steichiometry_dict]
    chains_list = list(itertools.chain.from_iterable(chains_list))

    #REMOVED THE RANDOM START TO A AND B. MAYBE THERE WOULD BE SOME FUTURE IMPLEMENTATION
    #random_start = random.sample(chains_list, 2)
    #random_start_1 = random.choice(list(pairwise_dict.keys()))
    #random_start_2 = random.choice(list(pairwise_dict[random_start_1].keys()))
    
    random_start_1 = "A"
    random_start_2 = "B"

    #if args['verbose']:
    #    print ("Starting from random chains:", random_start[0] ,"and", random_start[1])

    start_complex = pairwise_dict[random_start_1][random_start_2]
    

    pdb_complex = os.path.join(args['output_folder'], "complex.pdb")
    copyfile(start_complex, pdb_complex)

    chains_list.remove(random_start_1)
    chains_list.remove(random_start_2)

    while len(chains_list) > 0:
        try_chain = random.choice(chains_list)

        file_to_try = []

        while len(file_to_try) == 0:

            try_chain = random.choice(chains_list)
            file_to_try = check_chain_addition(pdb_complex , try_chain , pairwise_dict , steichiometry_dict)

        print("Trying to add" , try_chain)

        for file_pdb in file_to_try:
            chain_addition = build_complex(pdb_complex , file_pdb)

            if chain_addition:
                chains_list.remove(try_chain)
                print("Chain added successfully. Number of remaining chains:" , len(chains_list))
                if len(chains_list) == 0:
                    print("All chains have been added.")
                if len(chains_list) != 0:
                    print("The remaining chains to add are" , chains_list)
                break
    
    return os.path.join(args['output_folder'], "superimposition_complex.pdb")

