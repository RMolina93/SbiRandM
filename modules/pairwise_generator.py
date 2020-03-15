import os, sys, warnings
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import glob

warnings.filterwarnings("ignore")

def write_pairwise(pdb_file, pairwise_interactions, tmp):
    """
    input: PDB file, list of pairwise chain id that makes interactions, obtained with the function pairwise_generator
    output: This function writes the PairWise interactions of chains in different PDB files in the same folder as input.
    """

    split_call = f'./modules/scripts/PDBtoSplitChain2.pl -i {pdb_file} -o {tmp}/temp_{os.path.basename(pdb_file)[:-4]}_' 
    os.system(split_call)
    for chainA, chainB in pairwise_interactions:
       ## !!! WARNING AÑADIR EL PDBTOADDCHAIN PARA PONER CADENAS A Y B EN LA PROTEINA 
        os.system(f'./modules/scripts/PDBtoAddChain.pl -i {tmp}/temp_{os.path.basename(pdb_file)[:-4]}_{chainA}.pdb -o {tmp}/temp_{os.path.basename(pdb_file)[:-4]}_{chainA}_mod -c "A"')
        os.system(f'./modules/scripts/PDBtoAddChain.pl -i {tmp}/temp_{os.path.basename(pdb_file)[:-4]}_{chainB}.pdb -o {tmp}/temp_{os.path.basename(pdb_file)[:-4]}_{chainB}_mod -c "B"')
        os.system(f'cat {tmp}/temp_{os.path.basename(pdb_file)[:-4]}_{chainA}_mod.pdb > {tmp}/complex_{chainA}_{chainB}.pdb')
        os.system(f'echo "TER" >> {tmp}/complex_{chainA}_{chainB}.pdb')
        os.system(f'cat {tmp}/temp_{os.path.basename(pdb_file)[:-4]}_{chainB}_mod.pdb >> {tmp}/complex_{chainA}_{chainB}.pdb')
    filelist = glob.glob(f"{tmp}/temp_*.*")
    #for f in filelist:
    #    os.remove(f)


def check_chain_interaction(chainA, chainB, threshold):
    """
    This chain takes two BioPython PDB Chain objects, and return True is any atom CB of the first chain
    is closer than threshold to one atom CB of the second chain. Otherwise, returns False.
    """
    if chainA == chainB:
        #print ("Chains are equal, skipping")
        return False
    else:
        for residue_A in chainA:
            for residue_B in chainB:
                for atom_A in residue_A:
                    for atom_B in residue_B:
                        if atom_A.id == "CB" and atom_B.id == "CB":
                            if atom_A - atom_B < threshold: 
                                print ("Chains interact")
                                return True
                            else: continue
        return False
       

def pairwise_generator(pdb_file, tmp):
    """
    This function takes a PDB file, and decompose it into the chains objects of Biopython.
    Then it calls the function check_chain_interactions, and returns a list of tupples containing the ID of the chains that interact.
    """

    io = PDBIO()
    parser = PDBParser(PERMISSIVE=1)
    interacting_chains = list()
    structure = parser.get_structure('Complex', pdb_file)
    for model in structure:
        if len(model) == 1:
            return ("PDB only has one chain.")
        else:
            #print ("PDB has the next number of chains:", len(model))
            for chain in model:
                for second_chain in model:
                    if check_chain_interaction(chain, second_chain, 20): interacting_chains.append((chain.id, second_chain.id))
            
            unique_pairwise_interactions = ([v for k, v in enumerate(interacting_chains) if v[::-1] not in interacting_chains[:k]])
            write_pairwise(pdb_file, unique_pairwise_interactions, tmp)
                        
    return None


if __name__ == "__main__":
    pdb_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])))
    pdb_file = os.path.basename(sys.argv[1])
    pdb_full_path = os.path.join(pdb_path, pdb_file)
    try:
        TMP_folder = os.path.join(pdb_path, "TMP")
        os.mkdir(TMP_folder)
    except:
        print ("Directory is already created, skipping.")

    print (pdb_file)
    print(pairwise_generator(pdb_full_path, TMP_folder))

