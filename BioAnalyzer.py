'''Class to find potential coding sequences in segment of DNA. This class further determines the likelihood of a coding sequence near an 'ATG'/'AUG' 
(start codon) through consensus nucleotides nearby and upstream. Functions include: find_start_lst, find_stop_lst, find_codon_lst, find_cds_length, 
find_cds, find_cds_likelihood.'''
import FastaParser
from CODONS import CODONS

class BioAnalyzer():
    def __init__(self):
        pass
    
    #Finds all possible start sequences ('ATG') and returns their index as list.
    def find_start_lst(self):
        start = 'ATG'
        start_lst = []
        idx = 0
        while idx < len(self.DNA_sequence):
            idx = self.DNA_sequence.find(start, idx)
            if idx == -1:
                break
            else:
                start_lst.append(idx)
            idx += 3
        return start_lst
    
    #Finds all stop sequences and returns their index as list.
    def find_stop_lst(self):
        codons = self.codon_list()
        start = self.find_start_lst()
        stop_lst = []

        for i in range(len(start)):
            index = codons[i].index(CODONS['*'][0]) * 3 if CODONS['*'][0] in codons[i] else -1
            index1 = codons[i].index(CODONS['*'][1]) * 3 if CODONS['*'][1] in codons[i] else -1
            index2 = codons[i].index(CODONS['*'][2]) * 3 if CODONS['*'][2] in codons[i] else -1

            lst = sorted([index, index1, index2])
            
            for idx in lst:
                if idx != -1:
                    stop_lst.append(idx+start[i])
                    break
        return stop_lst

    #Finds a codon list based on start index.
    def codon_list(self):
        cdn_lst = []
        start = self.find_start_lst()
        end = len(self.RNA_sequence)
        for strt in start:
            lst = []
            for i in range(strt, end, 3):
                codon = self.RNA_sequence[i: i+3]
                if codon in CODONS['*']:
                    lst.append(codon)
                    break
                else:
                    lst.append(codon)
            cdn_lst.append(lst)

        return cdn_lst

    #Finds length of possible coding sequences and returns their length.
    def find_CDS_length(self):
        start = self.find_start_lst()
        stop = self.find_stop_lst()
        start = start[:len(stop)]
        CDS_length = [stop[i] - start[i] for i in range(len(start))]
        return CDS_length

    #Finds possible coding sequences and returns them as a list.
    def find_CDS(self):
        start = self.find_start_lst()
        stop = self.find_stop_lst()
        start = start[:len(stop)]
        CDS = [self.RNA_sequence[start[i]: stop[i]] for i in range(len(start))]
        return CDS
    
    #Algorithm to determine likelihood of coding sequence.
    def start_codon_likelihood(self):
        start = self.find_start_lst()
        stop = self.find_stop_lst()
        start = start[:len(stop)]
        start_codon_chance = []

        for strt in start:
            chance = 0
            if strt > 40 and 'TATA' in self.DNA_sequence[strt-40: strt-20]:
                chance += 1
            if strt > 6 and 'G' == self.DNA_sequence[strt-6] or 'C' == self.DNA_sequence[strt-6]:
                chance += 1
            if strt > 4 and 'G' == self.DNA_sequence[strt-4] or 'C' == self.DNA_sequence[strt-4] or 'A' == self.DNA_sequence[strt-4]:
                chance += 1
            if strt > 3 and 'G' == self.DNA_sequence[strt-3] or 'C' == self.DNA_sequence[strt-3] or 'A' == self.DNA_sequence[strt-3]:
                if self.DNA_sequence[strt-3] == 'A' or self.DNA_sequence[strt-3] == 'G':
                    chance += 2
                else:
                    chance += 1
            if strt > 2 and 'A' == self.DNA_sequence[strt-2] or 'C' == self.DNA_sequence[strt-2]:
                chance += 1
            if strt > 1 and 'G' == self.DNA_sequence[strt-1] or 'C' == self.DNA_sequence[strt-1] or 'A' == self.DNA_sequence[strt-1]:
                chance += 1
            if self.DNA_sequence[strt+3] == 'G' or self.DNA_sequence[strt+3] == 'A':
                chance += 1
            if self.DNA_sequence[strt+4] == 'C' or self.DNA_sequence[strt+4] == 'A':
                chance += 1
            start_codon_chance.append(chance)
        return start_codon_chance

    
if __name__ == '__main__':
    file1 = 'sequence.fasta'
    new_lst = fastaparser.parseFasta(file1)

    comp = []
    #Create test case
