'''Class to perform basic DNA to amino acid manipulations.'''
from random import choice
import fastaparser
from CODONS import CODONS

class BioTranslator:
    def __init__(self):
        pass
    
    #Return list of possible amino acid lengths.
    def amino_acid_length(self):
        nucleotides = self.find_CDS_length()
        AA_length = [DNA // 3 for DNA in nucleotides]
        return AA_length
        
    #Translate list of DNA coding sequences to list of amino acid sequences and return list.
    def translate(self):
        start = self.find_start_lst()
        stop = self.find_stop_lst()
        start = start[:len(stop)]
        codons_lst = list(CODONS.items())

        for j in range(len(start)):
            transl_seq = ''
            for i in range(start[j], stop[j], 3):
                codon = self.RNA_sequence[i: i+3]
                for aa in codons_lst:
                    if codon in aa[1]:
                        transl_seq += aa[0]
                        break
                    else:
                        continue
            self.AA_sequence.append(transl_seq)
        return self.AA_sequence

    #Reverse translate list of amino acids to list of RNA sequences.
    def rev_translate(self):
        nucleotides_lst = []

        for j in range(len(self.AA_sequence)):
            nucleotides = ''
            amino_acid_sequence = self.AA_sequence[j].upper()
            for i in range(0, len(amino_acid_sequence)):
                aa = amino_acid_sequence[i]
                #Choose which codon randomly from list of amino acid codons
                nucleotides += choice(CODONS[aa])
            nucleotides_lst.append(nucleotides)

        return nucleotides_lst
    
if __name__ == '__main__':
    #Create test case
    pass
    
