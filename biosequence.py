'''Basic class to perform bioinformatics analysis. Functions: complement, reverse complement, 
transcribe, translate, find start codon, find stop codons, find nucleotide length between start 
and stop, find amino acid length sequence alignments, gc content, and determine percent homology between 
DNA strands of the same or differing length'''
import fastaparser
from random import choice

CODONS = ({'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUA', 'CUU'], 'I': ['AUU', 'AUC', 'AUA'],
        'M': ['AUG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU',
        'AGC'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'],
        'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'Y': ['UAU', 'UAC'], '*': ['UAA', 'UAG',
        'UGA'], 'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAU', 'AAC'], 'K': 
        ['AAA', 'AAG'], 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['UGU', 'UGC'],
        'W': ['UGG'], 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG']})

class BioSequence:
    def __init__(self, DNA_sequence='', AA_sequence='', RNA_sequence='', gene='', organism=''):
        self.DNA_sequence = DNA_sequence.upper()
        self.AA_sequence = AA_sequence.upper()
        self.RNA_sequence = RNA_sequence.upper()
        self.DNA_complement_sequence = ''
        self.DNA_reverse_complement = ''
        self.get_DNA_RNA_AA()

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


    def find_start(self):
        start = 'ATG'
        return self.DNA_sequence.find(start)

    #Helper function for other functions needing to convert codons and amino acids.
    def codon_list(self):
        lst = []
        start = self.find_start()
        end = len(self.RNA_sequence)

        for i in range(start, end, 3):
            codon = self.RNA_sequence[i: i+3]
            if codon in CODONS['*']:
                lst.append(codon)
                break
            else:
                lst.append(codon)

        return lst

    #Helper function to return stop codon index.
    def find_stop(self):
        codons = self.codon_list()
        start = self.find_start()
        
        index = codons.index(CODONS['*'][0]) * 3 if CODONS['*'][0] in codons else -1
        index1 = codons.index(CODONS['*'][1]) * 3 if CODONS['*'][1] in codons else -1
        index2 = codons.index(CODONS['*'][2]) * 3 if CODONS['*'][2] in codons else -1

        lst = sorted([index, index1, index2])
        
        for idx in lst:
            if idx != -1:
                return idx + start
        return start

    def find_CDS_length(self):
        start = self.find_start()
        stop = self.find_stop()
        return int(stop - start)

    def find_CDS(self):
        start = self.find_start()
        stop = self.find_stop()
        return self.RNA_sequence[start: stop]

    def amino_acid_length(self):
        nucleotides = self.find_CDS_length()
        return nucleotides // 3

    def gc_content(self):
        gc = 0.0
        for nuc in self.DNA_sequence:
            if nuc == 'G' or nuc == 'C':
                gc += 1
        return gc / len(self.DNA_sequence)

    def translate(self):
        transl_seq = ''
        self.transcribe()
        start = self.find_start()
        stop = self.find_stop()

        for i in range(start, stop, 3):
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
        self.AA_sequence = transl_seq
        return transl_seq

    def rev_translate(self):
        nucleotides = ''
        amino_acid_sequence = self.AA_sequence.upper()

        for i in range(0, len(amino_acid_sequence)):
            aa = amino_acid_sequence[i]
            #Choose which codon randomly from list of amino acid codons
            nucleotides += choice(CODONS[aa])

        return nucleotides

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
            return
        if self.DNA_sequence:
            self.transcribe()
            self.translate()
            return

    def __repr__(self):
        DNA = self.DNA_sequence
        gc = str(self.gc_content())
        Comp = self.complement()
        Rcomp = self.rev_complement()
        RNA = self.RNA_sequence
        RTransc = self.rev_transcribe()
        CDS_length = str(self.find_CDS_length())
        CDS = self.find_CDS()
        AA = self.AA_sequence
        AA_length = self.amino_acid_length()
        RTransl = self.rev_translate()
        return ('''GC: %s DNA: %s \n\nComplement: %s \n\nReverse Complement: %s\nRNA: %s 
        \n\nReverse Transcribe: %s \n\nCDS Length: %s \nCoding Sequence:   %s \nReverse Translate: %s (variable tRNAs)
        \nAA Length: %s \nAmino Acid Sequence: %s\n''' % (gc, DNA, 
        Comp, Rcomp, RNA, RTransc, CDS_length, CDS, RTransl, AA_length, AA))


if __name__ == '__main__':
    file1 = 'sequence.fasta'
    new_lst = fastaparser.parseFasta(file1)

    comp = []
    for species in new_lst:
        seq = BioSequence(DNA_sequence=species[2])
        comp.append(seq)
    for sequence in comp:
        print(sequence)
    #print(comp[7].homology(comp[1].DNA_sequence))