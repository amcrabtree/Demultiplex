#!/usr/bin/env python3

import numpy as np
import re

# Constants:
DNA_bases = ("A", "T", "G", "C")
RNA_bases = ("A", "U", "G", "C")
test_seq = "TGCAGGTTGAGTTCTCGCTGTCCCGCCCCCATCTCTTCTCTTCCAGTCTGGCTCTGGAGCAGTTGAGCCCAGCTCAGGTCCGCCGGAGGAGACCG"
test_phred_score = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"

#####################################################################################

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter)-33

if __name__ == "__main__":
	assert convert_phred("A") == 32, "Phred score incorrect"
	assert convert_phred("@") == 31, "Phred score incorrect"
	assert convert_phred("#") == 2, "Phred score incorrect"
	print("PASS\tconvert_phred")

#####################################################################################

def qual_score(phred_score: str) -> float:
    '''Converts the phred score line from a read and calculates average quality score of read'''
    psum = 0
    plen = len(phred_score)
    for i in range(plen):
        psum += convert_phred(phred_score[i])
    av_q = psum/plen
    return av_q

if __name__ == "__main__":
	assert int(qual_score(test_phred_score)) == 37, "Average qual score incorrect"
	print("PASS\tqual_score")

#####################################################################################

def validate_base_seq(seq: str, RNAflag: bool = False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = set(seq.upper())
    model = set("AUGC" if RNAflag else "ATGC")
    return seq <= model  # if seq is within model, will be True

if __name__ == "__main__":
	assert validate_base_seq(
	    "AATAGAT") == True, "Validate base seq does not work on DNA"
	assert validate_base_seq(
	    "AAUAGAU", True) == True, "Validate base seq does not work on RNA"
	print("PASS\tvalidate_base_seq")

#####################################################################################

def gc_content(seq: str) -> float:
    '''Calculates the GC content of a DNA sequence'''
    seq = seq.upper()
    gc = seq.count("C") + seq.count("G")
    gc_cont = gc/len(seq)
    return gc_cont

if __name__ == "__main__":
	assert gc_content("GCGCGC") == 1
	assert gc_content("AATTATA") == 0
	assert gc_content("GCATGCAT") == 0.5
	print("PASS\tgc_content")

#####################################################################################

def oneline_fasta(infile: str, outfile: str):
    '''
    Converts input fasta file into nicer format fasta.
    Works by only adding newlines right before each header
    unless header is very first line in input file.
    '''
    ofh = open(outfile, "a")
    ofh.truncate(0)  # clears file before appending
    with open(infile, "r") as ifh:
        i = 0
        for line in ifh:
            i += 1
            line = line.strip("\n")
            if line.startswith(">"):
                if i != 1:
                    ofh.write("\n")
                ofh.write(line+"\n")
            else:
                ofh.write(line)
    ofh.close()
    print(infile, "converted to", outfile)

#####################################################################################
def read_stats(file: str) -> tuple:
    '''Returns average read length of input fastq file.'''
    with open(file, "r") as f:
        nt_count, max_rlen, min_rlen = 0, 0, 0
        i = 0
        for line in f:
            i += 1
            if i==2:
                min_rlen = len(line.strip()) # set a number for min len
            if i % 4 == 2: # if it's a read seq line
                read_len = len(line.strip())
                nt_count += read_len # add nucleotide count to counter
                if read_len > max_rlen:
                    max_rlen = read_len
                if read_len < min_rlen:
                    min_rlen = read_len
        readnum = i/4 # number of reads in file
        av_rlen = nt_count/readnum
        rstat_tup = (av_rlen, min_rlen, max_rlen)
        return rstat_tup

if __name__ == "__main__":
	assert read_stats("test.fastq") == (319, 214, 402)
	print("PASS\tread_stats")

#####################################################################################
