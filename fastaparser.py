'''Class to parse FASTA files to test my biotools class. This class will be updated to better parse FASTA files. 
This class is only used for testing purposes and will be later updated to more properly parse FASTA files.'''

def parseFasta(fastafile):
    description = ''
    sequence = ''
    lst = []
    desc_lst = []

    with open(fastafile, 'r') as nfile:
        for line in nfile:
            if line[0] == '>':
                lst.append((desc_lst, description, sequence))
                description = line.rstrip()
                desc_lst = description.split(' ')
                sequence = ''
            else:
                sequence += line.rstrip()
        lst.append((desc_lst, description, sequence))
    return lst[1:]
