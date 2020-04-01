import sys, os
from Bio import SeqIO
from Bio import pairwise2
import itertools, difflib
from sbiRandM.sbiRandM.data import Alphabet

def check_fasta_stoichometry(fasta_path):

    """
    @input - Fasta file of the Complex

    @output - Dictionary with the steichiometry of the complex.
              Name of the chain : 
                        Steichiometry - Absolute number that the chain appear in the Fasta.
                        Sequence - Sequence of the chain
    """

    steichiometry_dict = dict()
    chains = list(SeqIO.parse(fasta_path, "fasta"))

    index = 0

    for chain in chains:

        chain.id = Alphabet[index]
        steichiometry_dict[chain.id] = {"steichiometry" : 1, "sequence" : str(chain.seq)}
        index += 1

    for a,b in itertools.combinations(chains, 2):

        if difflib.SequenceMatcher(None, str(a.seq), str(b.seq)).ratio() > 0.70:
            #print ("Chains", a.id, b.id, "match!")
            try:
                steichiometry_dict[a.id]["steichiometry"] += 1
                del steichiometry_dict[b.id]
            except:
                continue
    
    return steichiometry_dict

if __name__ == "__main__":
    print( check_fasta_stoichometry(sys.argv[1]) )

