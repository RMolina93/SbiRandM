import argparse, warnings
from sbiRandM.sbiRandM import *



def mainSuperimp(args):
   if not os.path.isdir(args['output_folder']):
      if args.verbose:
         print (f'Output folder {args.output_folder} do not exist. Creating...')
      os.mkdir(args['output_folder'])
   if args.verbose:
      print ("Checking the steichiometry of the complex using the fatsa file.")
   steichiometry_dict = check_fasta_stoichometry(args['fasta_seq'])

   if args.verbose:
      print (f'Parsing the pairwise interactions structures from {args.folder}')
   pairwise_dict = obtain_pairwise_dict(steichiometry_dict , args['folder'])
   # print (pairwise_dict)
   if args.verbose:
      print (f'Starting to build the complex.')
   try:
      pdb_complex = execute_complex(steichiometry_dict , pairwise_dict , args)
   except ValueError:
      raise BadFastaException
   """

   while not model_validation(pdb_complex , args['fasta_seq']):
      print("Something went wrong, redoing model again.")
      try:
         pdb_complex = execute_complex(steichiometry_dict , pairwise_dict , args['output_folder'])
      except:
         print("Wow, something went VERY wrong! Redoing model")
         continue
         
         """
