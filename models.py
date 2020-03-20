import os, sys

class Chain():

    def __init__(self, name, sequence, whichChain=list()):
        self.name = name
        self.sequence = str(sequence)
        self.used = False
        self.whichChain=list()
    
    def used(self):
        self.used = True

    def dna_to_placeholder(self):
        """
        This function checks if chain is DNA or Protein, and converts nucleotides to placeholders.
        """
        if any(not ["T","G","A","C"] for c in str(self.sequence)):
            print ("Chain is Protein")
            return None
        else:
            new_seq = self.sequence.replace("T","Z")
            new_seq = new_seq.replace("G","X")
            new_seq = new_seq.replace("C","U")
            new_seq = new_seq.replace("A","O")
            self.sequence = new_seq
        return None


class Protein_Interaction():

    def __init__(self, name, path):
        self.name = name
        self.chains =list()
        self.path = path
        self.used = False

    def add_chain(self,chain):
        self.chains.append(chain)

    def used(self):
        self.used = True

class Query():

    def __init__(self, name):
        self.name = name
        self.chains =list()

    def add_chain(self,chain):
        self.chains.append(chain)

