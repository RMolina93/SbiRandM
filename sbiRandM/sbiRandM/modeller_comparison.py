import sys
sys.path.insert(0,"/usr/local/Cellar/modeller/9.23/modlib")
#from modeller import *
#from modeller.automodel import *
import glob, os, argparse, warnings
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO
from Bio import pairwise2
import shutil

from sbiRandM.sbiRandM.data import aminoacids
from sbiRandM.sbiRandM.models import Protein_Interaction, Chain, Query
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
    #print ("Len of record dict", len(record_dict))
    query = Query(name=(os.path.splitext(os.path.basename(fasta))[0]))
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    index = 0
    for sequence in record_dict.keys():
        print (record_dict[sequence], len(record_dict[sequence]))
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
    for pdb_file in glob.glob(os.path.join(folder,"*.pdb")):
        pdb_name = os.path.basename(pdb_file)
        structure = parser.get_structure('Complex', pdb_file)
        protein = Protein_Interaction(name = pdb_name, path = pdb_file)


        for pdb_chain in structure[0]:
            sequence = ""
            
            for residue in pdb_chain:
                sequence += aminoacids[residue.get_resname().strip()]
            protein.add_chain(Chain(name = pdb_chain.get_id(), sequence = sequence, first_aminoacid = list(pdb_chain)[0].get_id()[1], last_aminoacid = list(pdb_chain)[-1].get_id()[1]))

        list_of_interactions.append(protein)
    #print (list_of_interactions)
    return list_of_interactions

def check_similarity(query, interactions_list):

    """
    This function takes a query object and a list of Pairwise interaction objects,
    Returns the interaction list with the attribute "originalChain" updated for each chain of
    the pairwise interactions.  
    First it tries to check if the sequence of interaction is a subset. If not, it checks over a 73% of similarity
    """

    for chain in query.chains:

        for protein in interactions_list:
            for protein_chain in protein.chains:

                #Check if query is a subset
                if protein_chain.sequence in chain.sequence:  
                    protein_chain.originalChain.append(chain.name)

                # If not subset, check if its same chain with some discordance.
                elif difflib.SequenceMatcher(None,chain.sequence.strip(), str(protein_chain.sequence)).ratio() >0.73 :
                    protein_chain.originalChain.append(chain.name)

    return interactions_list

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
                    for letter in "ZXUOJ":
                        sequence = sequence.replace(letter,".")
                    output.write(sequence + "\n/\n")

                elif str(interaction.chains[1].originalChain[0]) == str(chain.name): # CHAIN B OF INTERACTION

                    alignments = pairwise2.align.globalms(chain.sequence, interaction.chains[1].sequence,  2, -1, -30, -10)
                    sequence = alignments[0][1]
                    for letter in "ZXUOJ":
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
            for letter in "ZXUOJ":
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

    a = automodel(env,
                alnfile  = pir_file, # alignment filename
                knowns   = tuple(templates),     # codes of the templates
                sequence = str(os.path.splitext(os.path.basename(fasta))[0].strip()))       

    a.starting_model= 1                 # index of the first model
    a.ending_model  = 2                 # index of the last model

    a.make()                            # do the actual homology modeling

def separe_interactions(interactions):
    """
    This function takes a list of protein pairwise interactions, that have more than originalChain
    (When mapping to homodimers) and divide them into several objects, with one chain for each.
    """
    updated_interactions = list()

    for interaction in interactions:
        for original_chain_A in interaction.chains[0].originalChain:
            for original_chain_B in interaction.chains[1].originalChain:
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
    """
    This function removes all the temporal files generated
    for the usage of Modeller
    """

    to_remove = glob.glob(os.path.join(args['folder'], "*SEPARED*"))
    to_remove = to_remove + glob.glob(os.path.join(args['output_folder'], "*SEPARED*")) + \
                            glob.glob(os.path.join(output_folder, "*.V999*")) + \
                            glob.glob(os.path.join(output_folder, "*.ini*")) + \
                            glob.glob(os.path.join(output_folder, "*.rsr*")) + \
                            glob.glob(os.path.join(output_folder, "*.D000*")) + \
                            glob.glob(os.path.join(output_folder, "*.sch*")) + \
                            glob.glob(os.path.join(output_folder, "*.pir*"))

    for file in to_remove:
        os.remove(file)



def mainMod(args):
    query = fasta_to_object(args['fasta_seq'])
    interactions = create_models(args['folder'])
    check_similarity(query, interactions)
    updated_interactions = separe_interactions(interactions)
    reorder_pdb(updated_interactions)
    
    generate_alignment(query, updated_interactions, args['output_folder'])
    
    make_model(args['output_folder'], args['folder'], args['fasta_seq'])
    clean_directories(args['output_folder'], args['folder'])
    



