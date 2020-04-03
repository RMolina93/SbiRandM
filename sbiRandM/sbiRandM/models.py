import os, sys

class Chain():
    """
    This object stores the properties of a polypeptide chain, 
    such as the name, the aminoacid sequence, the chain of the fasta file which is
    mapped to, and the first and last aminoacid numbers.
    """

    def __init__(self, name, sequence, first_aminoacid, last_aminoacid):
        self.name = name
        self.sequence = str(sequence)
        self.originalChain = list()
        self.first_aminoacid = str(first_aminoacid)
        self.last_aminoacid = str(last_aminoacid)


    def dna_to_placeholder(self):
        """
        This function checks if chain is DNA or Protein, and converts nucleotides to placeholders.
        """
        if all(c in "TGAC" for c in str(self.sequence)):
            new_seq = self.sequence.replace("T","Z")
            new_seq = new_seq.replace("G","X")
            new_seq = new_seq.replace("C","U")
            new_seq = new_seq.replace("A","O")
            self.sequence = new_seq
            return True
        else:
            return False

    def pretty_print(self):
        """
        This function prints the information of the Chain object in a more agradable way
        """
        print ("Chain name:", self.name)
        print ("Chain sequence:", self.sequence)
        print ("Chain Original Chain:" , self.originalChain)
        print ("First Aminoacid Number:", self.first_aminoacid)
        print ("Last Aminoacid Number:", self.last_aminoacid)

class Protein_Interaction():
    """
    This class defines an object of a Pairwise Interaction (Protein-Protein, Protein-DNA or DNA-DNA)
    It stores information such as the path of the structure, and if the chains are reversed in order
    according to the fasta.
    """

    def __init__(self, name, path):
        self.name = name
        self.chains =list()
        self.path = path
        self.reversed = False
        self.first_chain = "A"
        self.last_chain = "B"

    def add_chain(self,chain):
        self.chains.append(chain)

    def reverse(self):
        self.first_chain = "B"
        self.last_chain = "A"

    def pretty_print(self):
        print ("Interaction name:", self.name)
        print ("Chains:", self.chains)
        print ("Path:" , self.path)
        print ("Biopython object:", self.biopython_object)
        print ("Reversed:", self.reversed)


class Query():
    """
    This class is a container for the chains of the query fasta
    """

    def __init__(self, name):
        self.name = name
        self.chains =list()

    def add_chain(self,chain):
        self.chains.append(chain)

