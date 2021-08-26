'''Class to perform basic DNA to amino acid manipulations.'''
from CODONS import CODONS
from random import choice
import FastaParser

class BioTranslator:
    def __init__(self):
        pass
    
    #Return list of possible amino acid lengths.
    def amino_acid_length(self):
        nucleotides = self.find_CDS_length()
        AA_length = [DNA // 3 for DNA in nucleotides]
        return 
        
    #Translate list of DNA coding sequences to list of amino acid sequences and return list.
    def translate(self):
        start = self.find_start_lst()
        stop = self.find_stop_lst()
        start = start[:len(stop)]

        for j in range(len(start)):
            transl_seq = ''
            for i in range(start[j], stop[j], 3):
                codon = self.RNA_sequence[i: i+3]
                if start == -1:
                    break
                if codon in CODONS['F']:
                    transl_seq += 'F'
                elif codon in CODONS['L']:
                    transl_seq += 'L'
                elif codon in CODONS['I']:
                    transl_seq += 'I'
                elif codon in CODONS['M']:
                    transl_seq += 'M'
                elif codon in CODONS['V']:
                    transl_seq += 'V'
                elif codon in CODONS['S']:
                    transl_seq += 'S'
                elif codon in CODONS['P']:
                    transl_seq += 'P'
                elif codon in CODONS['T']:
                    transl_seq += 'T'
                elif codon in CODONS['A']:
                    transl_seq += 'A'
                elif codon in CODONS['Y']:
                    transl_seq += 'Y'
                elif codon in CODONS['*']:
                    transl_seq += '*'
                    break
                elif codon in CODONS['H']:
                    transl_seq += 'H'
                elif codon in CODONS['Q']:
                    transl_seq += 'Q'
                elif codon in CODONS['N']:
                    transl_seq += 'N'
                elif codon in CODONS['K']:
                    transl_seq += 'K'
                elif codon in CODONS['D']:
                    transl_seq += 'D'
                elif codon in CODONS['E']:
                    transl_seq += 'E'
                elif codon in CODONS['C']:
                    transl_seq += 'C'
                elif codon in CODONS['W']:
                    transl_seq += 'W'
                elif codon in CODONS['R']:
                    transl_seq += 'R'
                elif codon in CODONS['S']:
                    transl_seq += 'S'
                elif codon in CODONS['G']:
                    transl_seq += 'G'
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
    
