import sys
sys.path.insert(0,"/usr/local/Cellar/modeller/9.23/modlib")
from modeller import *
from modeller.automodel import *
import glob, os, argparse, warnings
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO
from Bio import pairwise2
import shutil

from .data import modeller_aminoacids as aminoacids
from .models import Protein_Interaction, Chain, Query
import difflib

warnings.filterwarnings("ignore")


def fasta_to_object(fasta):

    """
    This function takes a fasta file with several chains, and parse it into a
    query object. Also changes the DNA chains to Placeholder characters.

    @input folder - File with Fasta
    @output - Query object.
    """

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    #print (record_dict)
    #print ("Len of record dict", len(record_dict))
    query = Query(name=(os.path.splitext(os.path.basename(fasta))[0]))
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    index = 0
    for sequence in record_dict.keys():
        #print (record_dict[sequence], len(record_dict[sequence]))
        chain = Chain(name = alphabet[index], sequence = record_dict[sequence].seq, first_aminoacid = 1, last_aminoacid = len(record_dict[sequence].seq))
        chain.dna_to_placeholder()
        query.add_chain(chain)
        index += 1
    #print ("Fasta to object query:", query.chains)
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
    set_residues = set()

    for pdb_file in glob.glob(os.path.join(folder,"*.pdb")):
        pdb_name = os.path.basename(pdb_file)
        structure = parser.get_structure('Complex', pdb_file)
        protein = Protein_Interaction(name = pdb_name, path = pdb_file)


        for pdb_chain in structure[0]:
            sequence = ""
            for residue in pdb_chain:
                set_residues.add(residue.get_resname())
                sequence += aminoacids[residue.get_resname().strip()]
            protein.add_chain(Chain(name = pdb_chain.get_id(), sequence = sequence, first_aminoacid = list(pdb_chain)[0].get_id()[1], last_aminoacid = list(pdb_chain)[-1].get_id()[1]))

        list_of_interactions.append(protein)
        #print (protein.name)
    #print (set_residues)
    #for element in set_residues:
    #    print (element, aminoacids[element.strip()])

    #print (list_of_interactions)
    return list_of_interactions

def check_similarity(query, interactions_list):

    """
    This function takes a query object and a list of Pairwise interaction objects,
    Returns the interaction list with the attribute "originalChain" updated for each chain of
    the pairwise interactions.  
    First it tries to check if the sequence of interaction is a subset. If not, it checks over a 73% of similarity
    """
    print (len(query.chains))
    for chain in query.chains:
        #print (chain.name)
        for protein in interactions_list:

            for protein_chain in protein.chains:

                if protein_chain.sequence in chain.sequence:
                    if any(c in chain.sequence for c in "ZJ_OB"):
                        if len(chain.sequence) > len(protein_chain.sequence):
                            chain.sequence = protein_chain.sequence  
                    protein_chain.originalChain.append(chain.name)

                elif difflib.SequenceMatcher(None,chain.sequence.strip(), str(protein_chain.sequence)).ratio() >0.50 :
                    if any(c in chain.sequence for c in "ZJ_OB"):
                        if len(chain.sequence) > len(protein_chain.sequence):
                            chain.sequence = protein_chain.sequence 
                    protein_chain.originalChain.append(chain.name)
                


    return query, interactions_list

def generate_alignment(query, interactions, output_folder):

    if not os.path.isdir(output_folder): os.mkdir(output_folder)

    with open(os.path.join(output_folder, "alignment.pir"), "w") as output:

        # FIRST THE INTERACTIONS
        for interaction in interactions:
            output.write(">P1;" + interaction.name + "\n")
            output.write("structureX:" + 
                        interaction.name + ":" + 
                        interaction.chains[0].first_aminoacid + ":" + 
                        interaction.first_chain + ":" + 
                        interaction.chains[1].last_aminoacid + 
                        ":" + interaction.last_chain + ": : :-1.0:-1.0\n")
            
            for chain in query.chains:
                
                if str(interaction.chains[0].originalChain[0]) == str(chain.name): # CHAIN A OF INTERACTION

                    alignments = pairwise2.align.globalms(chain.sequence, interaction.chains[0].sequence,  2, -1, -30, -10)
                    sequence = alignments[0][1]
                    for letter in "ZB_OJ":
                        sequence = sequence.replace(letter,".")
                    output.write(sequence + "\n/\n")

                elif str(interaction.chains[1].originalChain[0]) == str(chain.name): # CHAIN B OF INTERACTION

                    alignments = pairwise2.align.globalms(chain.sequence, interaction.chains[1].sequence,  2, -1, -30, -10)
                    sequence = alignments[0][1]
                    for letter in "ZB_OJ":
                        sequence = sequence.replace(letter,".")
                    if chain.name == query.chains[-1].name:  
                        output.write(sequence + "*\n\n")
                    else:                    
                        output.write(sequence + "\n/\n")
                
                else:
                    if chain.name == query.chains[-1].name: 
                        output.write("-" * len(chain.sequence) + "*\n\n")
                    else:
                        output.write("-" * len(chain.sequence) + "\n/\n")


        #QUERY
        output.write(">P1;" + query.name.strip() + "\n")
        output.write("sequence:" + query.name + 
                     ": 1:" + query.chains[0].name + 
                     " :" + str(len(query.chains[0].sequence)) + 
                     ":" + query.chains[-1].name + 
                     ": : :-1.0:-1.0\n")
        
        for query_chain in query.chains:
            for letter in "ZB_OJ":
                query_chain.sequence = query_chain.sequence.replace(letter,".")

            if query_chain.name == query.chains[-1].name:  
                output.write(query_chain.sequence.strip() + "*\n\n") 
            else:                    
                output.write(query_chain.sequence.strip() + "\n/\n") 
            
def make_model(output_folder, interaction_pdb_folder, fasta):
    #detect pir file
    templates = list()
    for interaction_file in os.listdir(interaction_pdb_folder):
        if "SEPARED" in os.path.basename(interaction_file):
            templates.append(interaction_file)
            shutil.copyfile( os.path.join(interaction_pdb_folder, interaction_file) , os.path.join(output_folder, interaction_file) )

    os.chdir(output_folder)
    log.verbose()    # request verbose output
    env = environ()  # create a new MODELLER environment to build this model in

    # directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']
    pir_file = glob.glob("*.pir")[0]
    env.io.hetatm=True

    a = automodel(env,
                alnfile  = pir_file, # alignment filename
                knowns   = tuple(templates),     # codes of the templates
                sequence = str(os.path.splitext(os.path.basename(fasta))[0].strip()))       

    a.starting_model= 1                 # index of the first model
    a.ending_model  = 1                 # index of the last model

    a.make()                            # do the actual homology modeling

def separe_interactions(interactions):
    """
    This function takes a list of protein pairwise interactions, that have more than originalChain
    (When mapping to homodimers) and divide them into several objects, with one chain for each.
    """
    updated_interactions = list()

    for interaction in interactions:
        
        for original_chain_A in interaction.chains[0].originalChain:
            #print ("This complex enters the first for, original chain", interaction.chains[0].originalChain )
            for original_chain_B in interaction.chains[1].originalChain:
                print (interaction.name, original_chain_A, original_chain_B )
                #print ("This complex enters the second for, original chains", interaction.chains[1].originalChain )

                if original_chain_A == original_chain_B:
                    continue


                Updated_Chain_A = Chain(name = interaction.chains[0].name , 
                                        sequence = interaction.chains[0].sequence, 
                                        first_aminoacid = interaction.chains[0].first_aminoacid,
                                        last_aminoacid = interaction.chains[0].last_aminoacid)
                Updated_Chain_A.originalChain.append(original_chain_A)

                Updated_Chain_B = Chain(name = interaction.chains[1].name , 
                                        sequence = interaction.chains[1].sequence, 
                                        first_aminoacid = interaction.chains[1].first_aminoacid,
                                        last_aminoacid = interaction.chains[1].last_aminoacid)
                Updated_Chain_B.originalChain.append(original_chain_B)


                #COPY THE RESIDUE

                new_interaction_name = "SEPARED_" + interaction.name + "_" + original_chain_A + "_" + original_chain_B + ".pdb"
                new_path =  os.path.join(os.path.dirname(interaction.path), new_interaction_name)
                shutil.copyfile(interaction.path, new_path)

                Updated_interaction = Protein_Interaction(name = new_interaction_name,
                                                          path = new_path)

                Updated_interaction.add_chain(Updated_Chain_A)
                Updated_interaction.add_chain(Updated_Chain_B)

                if Updated_Chain_A.originalChain > Updated_Chain_B.originalChain:
                    Updated_interaction.reversed = True

                updated_interactions.append(Updated_interaction)

    return updated_interactions

def remove_duplicate_chains(folder):
    num_files = 0
    num_deleted_interactions = 0
    actual_list = list()
    for interaction_file in os.listdir(folder):
        num_files += 1
        if "SEPARED" in os.path.basename(interaction_file):
            if (interaction_file[-7], interaction_file[-5]) in actual_list or (interaction_file[-5], interaction_file[-7]) in actual_list:
                os.remove(os.path.join(folder,interaction_file))
                #print ("Deleted interaction:", interaction_file)
                num_deleted_interactions += 1
            else:
                actual_list.append((interaction_file[-7], interaction_file[-5])) 
    print ("Removed", num_deleted_interactions, "files from", num_files, ".")



def reorder_pdb(interactions):
    """
    This function takes a list of object interactions with the attribute Reversed, 
    and reorder the PDB files that has that attribute in True
    """

    for element in interactions:
        if element.reversed == True:
            element.reverse()
            element.chains = list(reversed(element.chains))

            corrected_PDB = os.path.join(os.path.dirname(element.path), "CORRECTED_" + element.name)

            with open(corrected_PDB, "w") as pdb_output:
                with open(element.path, "r") as pdb_input:
                    chain_A = list()

                    for line in pdb_input.readlines():
                        if "ATOM" in line and line[21] == "A":
                            chain_A.append(line)
                        elif "ATOM" in line and line[21] == "B":
                            pdb_output.write(line)
                    
                    pdb_output.write("TER\n")
                    
                    for line in chain_A:
                        pdb_output.write(line)

            os.remove(element.path)
            element.path = corrected_PDB
            element.name = os.path.basename(corrected_PDB)
            element.last_aminoacid = chain_A[-1].split(" ")[5]

            print ("Reversed PDB at path:", corrected_PDB)

def clean_directories(output_folder, TMP_folder):

    to_remove = glob.glob(os.path.join(TMP_folder, "*SEPARED*"))
    to_remove = to_remove + glob.glob(os.path.join(output_folder, "*SEPARED*")) + \
                            glob.glob(os.path.join(output_folder, "*.V999*")) + \
                            glob.glob(os.path.join(output_folder, "*.ini*")) + \
                            glob.glob(os.path.join(output_folder, "*.rsr*")) + \
                            glob.glob(os.path.join(output_folder, "*.D000*")) + \
                            glob.glob(os.path.join(output_folder, "*.sch*")) 

    for file in to_remove:
        os.remove(file)


def mainMod(args):

    query = fasta_to_object(args['fasta_seq'])
    interactions = create_models(args['folder'])
    query, interactions = check_similarity(query, interactions)
    updated_interactions = separe_interactions(interactions)
    reorder_pdb(updated_interactions)
    remove_duplicate_chains(args['folder'])
    
    generate_alignment(query, updated_interactions, args['output_folder'])
    
    make_model(args['output_folder'], args['folder'], args['fasta_seq'])

    #for interaction in interactions:
    #    interaction.pretty_print()

    clean_directories(args['output_folder'], args['folder'])
