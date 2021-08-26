# Biotools
IP: Need to create tests for each individual file. Tests are created for the BioTools.py file already.
    Will also update the FastaParser for more advanced parsing; for now, it is simply used for testing purposes.

BioTools uses composition and has a BioAnalyzer class object, BioSequence class object, and Biotranslator class object, and performs the following: Sequence manipulation, FASTA parsing, DNA homology, and determining possible coding sequences and their likelihood by examining consensus nucleotides.

There are eight files:
BioTools.py - main - Composition of BioAnalyzer, BioSequence, and BioTranslator
BioAnalyzer.py
BioSequence.py
BioTranslator.py
CODONS.py
FastaParser.py
sequence.fasta
README.md

To test the code, run the BioTools file, which has test code at the bottom. The Fastaparser is simple and designed only for testing purposes; it will be refined at a later time.

