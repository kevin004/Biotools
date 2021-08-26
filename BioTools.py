'''Main class that can handle DNA/RNA/AA manipulation, determine potential coding sequences and 
their likelihood based on consensus sequences, and analyze homology. Derived from BioSequence, BioTranslator, and BioAnalyzer.'''
from BioSequence import BioSequence
from BioTranslator import BioTranslator
from BioAnalyzer import BioAnalyzer
from CODONS import CODONS
import FastaParser

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
        if self.AA_sequence:
            self.rev_translate()
            self.rev_transcribe()
            return
        if self.DNA_sequence:
            self.transcribe()
            self.translate()
            return

    #Test each function of the classes through printing. 
    #Need to make this easier to read and more user friendly.
    def __repr__(self):
        DNA = self.DNA_sequence
        gc = str(self.gc_content())
        Comp = self.complement()
        Rcomp = self.rev_complement()
        RNA = self.RNA_sequence
        RTransc = self.rev_transcribe()
        CDS_length = str(self.find_CDS_length())
        CDS = self.find_CDS()
        SCL = self.start_codon_likelihood()
        AA = self.AA_sequence
        AA_length = self.amino_acid_length()
        RTransl = self.rev_translate()
        return ('''GC: %s DNA: %s \n\nComplement: %s \n\nReverse Complement: %s\nRNA: %s 
        \n\nReverse Transcribe: %s \n\nCDS Length: %s \nCoding Sequence:   %s SCL: %s
        \nReverse Translate: %s (variable tRNAs)
        \nAA Length: %s \nAmino Acid Sequence: %s\n''' % (gc, DNA, 
        Comp, Rcomp, RNA, RTransc, CDS_length, CDS, SCL, RTransl, AA_length, AA))
    


if __name__ == '__main__':
    file1 = 'sequence.fasta'
    new_lst = FastaParser.parseFasta(file1)

    comp = []
    for species in new_lst:
        seq = BioTools(DNA_sequence=species[2])
        comp.append(seq)
    for DNA in comp:
        print(DNA)
