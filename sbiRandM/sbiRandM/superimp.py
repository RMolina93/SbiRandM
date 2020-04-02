import argparse, warnings
from sbiRandM.sbiRandM import *



def mainSuperimp(args):
   if not os.path.isdir(args['output_folder']):
      os.mkdir(args['output_folder'])
   steichiometry_dict = check_fasta_stoichometry(args['fasta_seq'])
   pairwise_dict = obtain_pairwise_dict(steichiometry_dict , args['folder'])
   # print (pairwise_dict)
   pdb_complex = execute_complex(steichiometry_dict , pairwise_dict , args['output_folder'])
   
   while not model_validation(pdb_complex , args['fasta_seq']):
      print("Something went wrong, redoing model again.")
      try:
         pdb_complex = execute_complex(steichiometry_dict , pairwise_dict , args['output_folder'])
      except:
         print("Wow, something went VERY wrong! Redoing model")
         continue
