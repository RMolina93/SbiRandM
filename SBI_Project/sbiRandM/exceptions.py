
class FilesDontMatchException(Exception):

   def __init__(self): 
      pass
   def __str__(self):
      print ("The program has stopped because the files do not match the requirements.")
      return "The files do not match the requirements"

class FastaRaroException(Exception):

   def __init__(self): 
      pass
   def __str__(self):
      print ("The program stopped because the fasta file does not match with the files in the provided PDB folder")
      return "The fasta file does not match with the pdb directory provided"
   
class BadFastaException(Exception):

   def __init__(self): 
      pass
   def __str__(self):
      print ("The program stopped because your Fasta file is either bad configured or have loosen information.")
      return "Your fasta file is either bad configured or it's likely to have loosen information"