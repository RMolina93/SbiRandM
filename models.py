import os, sys

class Chain():

    def __init__(self, name, sequence, first_residue, last_residue, originalChain=list()):
        self.name = name
        self.sequence = str(sequence)
        self.used = False
        self.originalChain=list()
        self.first_aminoacid = str(first_residue)
        self.last_aminoacid = str(last_residue)

    def used(self):
        self.used = True

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

