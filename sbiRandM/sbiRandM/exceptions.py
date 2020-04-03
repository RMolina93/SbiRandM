
class FilesDontMatchException(Exception):
   "The files do not match the requirements"
   pass

class FastaRaroException(Exception):
   "The fasta file does not match with the pdb directory provided"
   pass

class BadFastaException(Exception):
   "Your fasta file is either bad configured or it's likely to have loosen information"
   pass