#!/usr/bin/env python3

import argparse

# define command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Normalizes the k-mer coverage of your input reads and saves those reads to a new FASTQ file")
	parser.add_argument("-1", "--R1_file", help="R1 (forward read) file", required=True)
	parser.add_argument("-2", "--R2_file", help="R2 (forward index) file", required=True)
	parser.add_argument("-3", "--R3_file", help="R3 (reverse index) file", required=True) 
	parser.add_argument("-4", "--R4_file", help="R4 (reverse read) file", required=True)
	parser.add_argument("-b", "--barcode_file", help="file of barcodes", required=True)
	parser.add_argument("-q", "--min_qscore", help="minimum qscore allowed", required=False, type=int)
	parser.add_argument("-o", "--output_dir", help="output folder for fastq files", required=True)
	return parser.parse_args()

# Assign variables from arguments in the command line
args = get_args()
kmer_size = args.kmer_len
f_reads = args.R1_file
r_reads = args.R4_file
f_index = args.R2_file
r_index = args.R3_file
dir_out = args.output_dir

# in this dictionary, kmers are keys and value is number of occurences of the kmer
fa_dict = {}

fho = open(fout, "a") # open output file
i: int=0 # keeps track of line number in fasta file
with open(filename, "r") as fh:
    reclist=[] 
    for line in fh:
        i+=1
        if i%4 == 1:
            reclist=[] # clear list of record info
            line=line.strip("\n")
            reclist.append(line)
        elif i%4 == 2:
            line=line.strip("\n")
            reclist.append(line)
        elif i%4 == 3:
            line=line.strip("\n")
            reclist.append(line)
        elif i%4 == 0:
            line=line.strip("\n")
            reclist.append(line)
            klst=[] # initialize list to store kmer cov scores
            # add kmer counts to main dictionary
            for n in range(len(reclist[1])-kmer_size+1):
                k=reclist[1][n:n+kmer_size]
                if k in fa_dict: