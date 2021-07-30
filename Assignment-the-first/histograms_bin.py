#!/usr/bin/env python

import numpy as np
import Bioinfo
import argparse
import gzip

# define command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Plots average qscore for all reads in input FASTQ file")
	parser.add_argument("-f", "--fq_in", help="input fastq file", required=True)
	parser.add_argument("-r", "--read_len", help="read length", required=True, type=int)
	parser.add_argument("-o", "--output", help="output file prefix", required=True)
	return parser.parse_args()

# Assign variables from arguments in the command line
args = get_args()
fq_in = args.fq_in
rlen = args.read_len
out_pic = args.output

## Use this cell to populate your numpy array 
#numreads=int(sum(1 for line in gzip.open(fq_in, "r"))/4)
def mean_qscores(fq_in, readlen) -> np.ndarray:
    av_npa = np.zeros((readlen), dtype=int)
    with gzip.open(fq_in,"r") as fh:
        i, r = 0, 0
        for line in fh:
            i+=1
            if i%4 == 0:
                r+=1
                line=line.rstrip()
                for basepos in range(len(line)):
                    p=line[basepos]-33
                    av_npa[basepos]+=p
        av_npa=av_npa/r
    return av_npa

mean=mean_qscores(fq_in, rlen)

print("position\tmean")
for i in range(len(mean)):
    print(i+1,mean[i],sep="\t")

## print plot
import matplotlib.pyplot as plt
xval=[i for i in range(1, rlen+1)]
plt.bar(xval, mean, color='blue')
plt.title("QScore Mean")
plt.xlabel("Base Pair Position")
plt.ylabel("Mean Quality Score")
plt.savefig(out_pic, format="png")