import argparse, warnings
from sbiRandM.sbiRandM import *

def mainSuperimp(args):
   if args['verbose']:
      print ("Verbose Mode: Activated.")

   if not os.path.isdir(args['output_folder']):
      if args['verbose']:
         print ('Output folder do not exist. Creating...')
      os.mkdir(args['output_folder'])

   if args['verbose']:
      print ("Checking the steichiometry of the complex using the fasta file.")
   steichiometry_dict = check_fasta_stoichometry(args['fasta_seq'])
   
   if args['verbose']:
      print ("The input fasta consist of", len(steichiometry_dict.keys()), "different chains.")

   if args['verbose']:
      print ('Parsing the pairwise interactions structures from', args['folder'])

   pairwise_dict = obtain_pairwise_dict(steichiometry_dict , args['folder'])
   if args['verbose']:
      print ("Parsed", len(pairwise_dict.keys()), "interaction structures.")

   if args['verbose']:
      print ('Starting to build the complex.')

   try:
      pdb_complex = execute_complex(steichiometry_dict , pairwise_dict , args)
      print ("The program has finished, and your complex is created in " + os.path.join(args['output_folder'],"superimposition_complex.pdb\n\n"))
   except ValueError:
      raise BadFastaException
