from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB import PDBParser
from Bio.PDB import Polypeptide
import numpy as np
import subprocess
import modeller_comparison as mc


pdbExtension = ".pdb"
proteinName = "1gzx"
proteinPath = "/Users/miguel/PycharmProjects/Structural_Project/structures/test/"
scriptPath = "/Users/miguel/PycharmProjects/Structural_Project/modules/scripts/PDBtoSplitChain2.pl"
stdoutPath = "/Users/miguel/PycharmProjects/Structural_Project/"


def superimpose():
   

   respuesta = mc.obtain_sequence('structures/6gmh/TMP/')
   
   str1 = list(respuesta.keys())[0]
   str2 = list(respuesta.keys())[0]

   sup = SVDSuperimposer()
   
   sup.set_atoms
if __name__ == "__main__":
   superimpose()
   
   
   '''
   if options.verbose:
      print('\nSuperimposing %s over chain %s' % (str2 , node.get_chain()))
   chain1 = node.get_chain()
   str2_copy = [x.copy() for x in str2]
   node_chain_copy = str1[0][chain1].copy()
   chain_str2 = i.copy()
   
   # Trim the chains so that they have the same number of atoms to superimpose
   trim_to_superimpose(node_chain_copy , chain_str2)
   atoms_chain1 = [atom for atom in list(node_chain_copy.get_atoms()) if atom.get_id() == 'CA' or atom.get_id() == 'P']
   atoms_chain2 = [atom for atom in list(chain_str2.get_atoms()) if atom.get_id() == 'CA' or atom.get_id() == 'P']
   sup = Superimposer()
   
   # Superimpose the chains
   sup.set_atoms(atoms_chain1 , atoms_chain2)
   
   # Select the chains that haven't been superimposed
   if not homodimer:
      other_chain2 = [x for x in str2_copy if x.get_id() != chain_str2.get_id()][0]
      other_chain2_original = [x for x in str2 if x.get_id() != chain_str2.get_id()][0]
   else:
      other_chain2 = str2_copy[1]
      other_chain2_original = str2[1]
   
   # Apply the rotation matrix to the chain we want to add to the complex
   sup.apply(other_chain2)
   
   # Assess if there is a clash between the chain we want to add and the others
   if not get_clash_chains(str1 , other_chain2 , chain1 , options):
      
      # add the chain to the macrocomplex
      other_chain2.id = len(complex_id.get_nodes()) + 1
      similar_seq[other_chain2] = similar_seq[other_chain2_original]
      complex_id.add_node(other_chain2 , node , str2)
      str1[0].add(other_chain2)
      
      if options.intensive:
         interaction_finder(str1 , other_chain2.get_id() , complex_id , node , options)
      
      return True
   else:
      return 'clash'
'''
'''
def trim_to_superimpose(chain1 , chain2):
   # Todo: mirar que pasa si consideras una similaridad del 100% para que sean similar seqs
   """
   Takes two chains and removes the residues that do not have a match in the sequence alignment.
   :param chain1: first chain.
   :param chain2:  second chain.
   :return: returns chain 1 and 2 trimmed so that they have the same atom length.
   """
   seq1 = get_sequence_from_chain(chain1)
   seq2 = get_sequence_from_chain(chain2)
   
   alignment = pairwise2.align.globalxx(seq1 , seq2)
   
   score = alignment[0][2]
   length = max(len(seq1) , len(seq2))
   ident_perc = score / length
   
   # in principle, there will not get here two chains that do not have more than 95% of similarity
   if ident_perc > 0.95:
      
      seq1_array = list(alignment[0][0])  # Storing the alignment sequences as arrays
      seq2_array = list(alignment[0][1])
      seq1_numeric = get_numeric_array(alignment[0][0])  # Storing those sequences with numbers instead of residuals
      seq2_numeric = get_numeric_array(alignment[0][1])
      to_delete_from_1 = []  # To avoid modifying a list while iterating through it we will store the elements to
      to_delete_from_2 = []  # remove in these lists
      pairs1 = zip(seq1_array , seq2_numeric)  # we pair the sequnce alignment with the numeric list of the other chain
      for pair in pairs1:
         if pair[0] == '-':  # If there is a gap in the sequence, we will remove the nth residue of the other chain
            to_delete_from_2.append(list(chain2.get_residues())[pair[1]].get_id())
      
      pairs2 = zip(seq2_array , seq1_numeric)
      for pair in pairs2:
         if pair[0] == '-':
            to_delete_from_1.append(list(chain1.get_residues())[pair[1]].get_id())
      
      for residue_to_delete in to_delete_from_1:  # Removing the residuals from the sequences
         chain1.__delitem__(residue_to_delete)
      
      for residue_to_delete in to_delete_from_2:  # Removing the residuals from the sequences
         chain2.__delitem__(residue_to_delete)
         '''