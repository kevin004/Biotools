'''Class to perform basic DNA/RNA manipulation. Functions: complement, reverse complement, 
transcribe,gc content.'''
import FastaParser
from CODONS import CODONS

class BioSequence:
    def __init__(self):
        pass

    def complement(self):
        comp_seq = ''

        for nucleotide in self.DNA_sequence:
            if nucleotide == 'A':
                comp_seq += 'T'
            elif nucleotide == 'T':
                comp_seq += 'A'
            elif nucleotide == 'G':
                comp_seq += 'C'
            elif nucleotide == 'C':
                comp_seq += 'G'
        self.DNA_complement_sequence = comp_seq
        return comp_seq

    def rev_complement(self):
        rev_comp = self.DNA_sequence[::-1]
        comp_seq = ''

        for nucleotide in rev_comp:
            if nucleotide == 'A':
                comp_seq += 'T'
            elif nucleotide == 'T':
                comp_seq += 'A'
            elif nucleotide == 'G':
                comp_seq += 'C'
            elif nucleotide == 'C':
                comp_seq += 'G'
        self.DNA_reverse_complement = comp_seq
        return comp_seq

    def transcribe(self):
        tran_seq = self.DNA_sequence.replace('T', 'U')
        self.RNA_sequence = tran_seq
        return tran_seq

    def rev_transcribe(self):
        rtran_seq = self.RNA_sequence.replace('U', 'T')
        self.DNA_sequence = rtran_seq
        return rtran_seq

    def gc_content(self):
        gc = 0.0
        for nuc in self.DNA_sequence:
            if nuc == 'G' or nuc == 'C':
                gc += 1
        return gc / len(self.DNA_sequence)

    def homology(self, sequence2):
        sequence1 = self.DNA_sequence
        homology = 0.0
        match = 0.0

        if len(sequence1) > len(sequence2):
            diff = len(sequence1) - len(sequence2)
            for i in range(diff):
                for j in range(i, len(sequence2) + i):
                    if sequence1[j] == sequence2[j-i]:
                        match += 1
                if match / len(sequence2) > homology:
                    homology = match / len(sequence2)
                match = 0.0
        elif len(sequence2) > len(sequence1):
            diff = len(sequence2) - len(sequence1)
            for i in range(diff):
                for j in range(i, len(sequence1) + i):
                    if sequence1[j-i] == sequence2[j]:
                        match += 1
                if match / len(sequence1) > homology:
                    homology = match / len(sequence1)
                match = 0.0
        else:
            for i in range(len(sequence1)):
                if sequence1[i] == sequence2[i]:
                    match += 1
            homology = match / len(sequence1)

        return homology


if __name__ == '__main__':
    file1 = 'sequence.fasta'
    new_lst = fastaparser.parseFasta(file1)

    comp = []
    for species in new_lst:
        seq = BioSequence(DNA_sequence=species[2])
        comp.append(seq)
    for DNA in comp:
        print(DNA)
