'''BioTools is a composition class that has a BioSequence, BioTranslator, and BioAnalyzer class object. It can handle DNA/RNA/AA manipulation, determine 
potential coding sequences and their likelihood based on consensus sequences, and analyze homology.'''
from biosequence import BioSequence
from biotranslator import BioTranslator
from bioanalyzer import BioAnalyzer
from CODONS import CODONS
import fastaparser

class BioTools(BioSequence, BioTranslator, BioAnalyzer):
    def __init__(self, DNA_sequence='', AA_sequence='', RNA_sequence='', gene='', organism=''):
        self.DNA_sequence = DNA_sequence.upper()
        self.AA_sequence = []
        if AA_sequence:
            self.AA_sequence += AA_sequence
        self.RNA_sequence = RNA_sequence.upper()
        self.DNA_complement_sequence = ''
        self.DNA_reverse_complement = ''
        self.get_DNA_RNA_AA()

    #Helper function to set initial variables.
    def get_DNA_RNA_AA(self):
        if self.DNA_sequence and self.RNA_sequence and self.AA_sequence:
            return
        if self.RNA_sequence:
            self.rev_transcribe()
            self.translate()
            return
        #Need to modify.
        if self.AA_sequence:
            self.rev_translate()
            return
        if self.DNA_sequence:
            self.transcribe()
            self.translate()
            return

    #Test each function of the classes through printing. 
    #Need to make this easier to read and more user friendly.
    def __repr__(self):
        DNA = self.DNA_sequence
        RNA = self.RNA_sequence
        RTransc = self.rev_transcribe()
        CDS = self.find_CDS()
        SCL = self.start_codon_likelihood()
        AA = self.AA_sequence
        RTransl = self.rev_translate()
        return ('''DNA: %s \nRNA: %s \nReverseTranscribed: %s
                \nCoding Sequences: %s \nStart Codon Likelihood: %s
                \nAmino Acid Sequence: %s \nReverse Translation: %s\n''' 
                % (DNA, RNA, RTransc, CDS, SCL, AA, RTransl))
    


if __name__ == '__main__':
    file1 = 'sequence.fasta'
    new_lst = fastaparser.parseFasta(file1)

    comp = []
    for species in new_lst:
        seq = BioTools(DNA_sequence=species[2])
        comp.append(seq)
    for DNA in comp:
        print(DNA)
