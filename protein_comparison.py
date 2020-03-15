from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB import PDBParser
from Bio.PDB import Polypeptide
import numpy as np
import subprocess


proteinPath = "/Users/miguel/PycharmProjects/Structural_Project/structures/test/1gzx.pdb"
scriptPath = "/Users/miguel/PycharmProjects/Structural_Project/modules/scripts/PDBtoSplitChain2.pl"

def superimposition():
   structure = PDBParser().get_structure("3t72", "structures/3t72.pdb")
   
   for x in structure.get_atoms():
      print(x.get_name())

   ref_model = structure[0]
   for alt_model in structure:
      ref_atoms = []
      alt_atoms = []
      for (ref_chain , alt_chain) in zip(ref_model , alt_model):
         for ref_res , alt_res , amino , allow in zip(ref_chain , alt_chain , seq_str , use):
            assert ref_res.resname == alt_res.resname
            assert ref_res.id == alt_res.id
            assert amino == Polypeptide.three_to_one(ref_res.resname)
            if allow:
               # CA = alpha carbon
               ref_atoms.append(ref_res['CA'])
               alt_atoms.append(alt_res['CA'])

      # Align these paired atom lists:
      super_imposer = SVDSuperimposer()
      super_imposer.set_atoms(ref_atoms , alt_atoms)
   
      if ref_model.id == alt_model.id:
         # Check for self/self get zero RMS, zero translation
         # and identity matrix for the rotation.
         assert np.abs(super_imposer.rms) < 0.0000001
         assert np.max(np.abs(super_imposer.rotran[1])) < 0.000001
         assert np.max(np.abs(super_imposer.rotran[0]) - np.identity(3)) < 0.000001
      else:
         # Update the structure by moving all the atoms in
         # this model (not just the ones used for the alignment)
         super_imposer.apply(alt_model.get_atoms())


#superimposition()
def decomposeInChains(proteinPath, scriptPath):
   #This function uses the perl script to decompose a protein in its chains
   
   # Result: Fastas and PDBs of the different chains
   # Format: <name_Of_Protein><Chain (A,B...)>.<fa>/<pdb>
   
   var = proteinPath
   pipe = subprocess.Popen(["perl", scriptPath ,"-i" ,var, "-o", "/Users/miguel/PycharmProjects/Structural_Project/1gzx"],stdout=subprocess.PIPE)
   pipe.communicate()