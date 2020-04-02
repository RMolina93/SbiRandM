import sys, os, warnings
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from Bio import SeqIO
from sbiRandM.sbiRandM import *
from sbiRandM.sbiRandM.data import Alphabet
import itertools, random
from shutil import copyfile


def Compute_equal_chain(structure_1, structure_2):
    """
    @ Input - Two Biopython PDB structures of an interaction.
    @ Output - A list of Atoms that correspond to the same chain.
    """

    for chain in structure_1[0]:
        for chain_2 in structure_2[0]:
            if Get_fasta(chain) == Get_fasta(chain_2):
                return (list(chain.get_atoms()), list(chain_2.get_atoms()))

    return (None, None)

def check_chain_addition(complex_pdb, chain, pairwise_dict, steichiometry_dict):
    """
    TO-DO
    Esto tiene que mirar si el complex PDB y la chain.id que queremos añadir, tienen alguna interaccion posible
    si la tienen, se devuelve el path del file
    si no la tiene, se devuelve false, por lo que se hace otro pick at random
    """
    parser = PDBParser(PERMISSIVE=1)

    complex_pdb = parser.get_structure('Complex', complex_pdb)
    
    chains_to_test_in_complex = list(pairwise_dict[chain].keys())

    common_chains = list()
    
    
    for chain_pdb in complex_pdb.get_chains():
        fasta = Get_fasta(chain_pdb)
        for chain_2 in chains_to_test_in_complex:
            if check_homology(fasta, steichiometry_dict[chain_2]["sequence"]):
                common_chains.append(pairwise_dict[chain][chain_2])
    
    return common_chains

    return None

def build_complex(file_1, file_2):
    # COMPLEX FILE, PAIRWISE_DICT AND CHAIN TO ADD???????
    """
    This function takes the complex output file (or in the first iteration one of the pairwise interactions)
    and another pairwise interaction PDB complex. Then it tries to add the chain to the complex until there is not clash

    @ Input - Two file path for a PDB interactions.
    @ Output - File path of the complex PDB file / Error: Chain cannot be added.
    """
    #print ("Trying to add", file_1, "and", file_2)

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
    #print(sup.rms)
    sup.apply(list(structure_2.get_atoms()))

    for chain in structure_2[0].get_chains():
        if chain.id != list(atoms_moving)[0].get_full_id()[2]:
            moved_chain = chain
            #print ("Going to move chain:", moved_chain.id)

    if check_clash(structure_1, moved_chain):
    #if True:
        with open(file_1, "wt") as out_file:
            
            for model in list(structure_1.get_chains()) + [moved_chain] :
                #print ("Hi im a model to write!")
                io.set_structure(model)
                io.save(out_file)

        rename_complex_chains(file_1)

        return True
    return False


    
    print ("\nAdded Chain", moved_chain.id, "to chains:", [i.id for i in structure_1.get_chains()])
    return output_name

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
    #print ("The PDB file is located in:", file)

def execute_complex(steichiometry_dict, pairwise_dict, output_folder):

    chains_list = [[i] * steichiometry_dict[i]["steichiometry"] for i  in steichiometry_dict]
    chains_list = list(itertools.chain.from_iterable(chains_list))

    random_start = random.sample(chains_list, 2)
    random_start_1 = random.choice(list(pairwise_dict.keys()))
    random_start_2 = random.choice(list(pairwise_dict[random_start_1].keys()))

    #print ("Starting from random chains:", random_start[0] ,"and", random_start[1])

    start_complex = pairwise_dict[random_start_1][random_start_2]
    pdb_complex = os.path.join(output_folder, "complex.pdb")
    copyfile(start_complex, pdb_complex)
    #print ("The PDB file is located in:", start_complex)

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
                print("The remaining chains to add are" , chains_list)

                break


        #print ("The remaining chains to add are", chains_list)
    
    return os.path.join(output_folder, "complex.pdb")

def model_validation(pdb_file, fasta):

    """
    Since this is a randomized algorithm, sometimes things could be wrong.
    Therefore we implemented a function, that redo the algorithm if the resulting PDB
    has a different number of chains than the fasta
    """

    parser = PDBParser(PERMISSIVE=1)
    complex_pdb = parser.get_structure('Complex', pdb_file)
    num_pdb_chains = len(list(complex_pdb.get_chains()))
    num_fasta_chains = len(list(SeqIO.parse(fasta, "fasta")))

    if num_fasta_chains == num_pdb_chains:
        return True
    else:
        return False



