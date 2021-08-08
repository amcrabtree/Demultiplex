#!/usr/bin/env python3

import argparse
import gzip
import Bioinfo

########################## ARGPARSE ##############################
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
f_reads = args.R1_file
r_reads = args.R4_file
f_index = args.R2_file
r_index = args.R3_file
barf = args.barcode_file
minq = args.min_qscore
dir_out = args.output_dir

########################## FUNCTIONS ##############################

def revc (seq: str) -> str:
    '''converts DNA seq to reverse complement'''
    seq=seq.upper()
    nucd = {"A":  "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    rseq=seq.translate(str.maketrans(nucd))[::-1] 
    return rseq

# def bc_correct(qbar_seq: str, barcode_txt: tuple):
#     '''input is barcode sequence in question and a tuple of possible
#     barcode values, output is corrected barcode or empty string 
#     if no correction possible'''
#     - determine min hamming dist in barcode tuple and save as variable
#         - initialize hamming dist variable as barcode seq length
#         - determine hamming dist between each barcode pair
#         - if dist < variable, replace
#     - open dust bin files (forward and reverse reads)

########################## MAIN ##############################

## open barcode file, store barcodes as dictionary, close file. 
## in bar_dict, value = index id, key = tuple of barcode & revcomp
b = open(barf, "r")
bar_dict={}
b.readline() # ignore first line
line = b.readline().rstrip()
while line:
    (sample,group,treatment,index,iseq) = line.split("\t")
    bar_dict[index] = iseq
    line = b.readline().rstrip()
b.close()

## open input files
fr = open(f_reads, "r")
rr = open(r_reads, "r")
fi = open(f_index, "r")
ri = open(r_index, "r")

'''
## open zipped input files
fr = gzip.open(f_reads, "rt")
rr = gzip.open(r_reads, "rt")
fi = gzip.open(f_index, "rt")
ri = gzip.open(r_index, "rt")
'''

## count lines in file (will be same for all files)
numrecs = int((max(i for i, _ in enumerate(ri))+1)/4)
ri = open(r_index, "r")

## generate/open output files (sample bins [# of barcodes *2], 
## swapping bins [2], dust bins [2])
swap_fw = open (dir_out+"/FW_swap", "a")
swap_rv = open (dir_out+"/RV_swap", "a")
junk_fw = open (dir_out+"/FW_junk", "a")
junk_rv = open (dir_out+"/RV_junk", "a")

## make dictionary of output file info while opening them
file_dict = {}
for el in bar_dict:
    f = dir_out + "/FW_" + el + ".fq"
    r = dir_out + "/RV_" + el + ".fq"
    fwh = "FW_"+el
    rvh = "RV_"+el
    file_dict[fwh] = open(f, "w")
    file_dict[rvh] = open(r, "w")

## record info loop
for i in range(numrecs):
    r1_lst = [next(fr).rstrip() for x in range(4)]
    r2_lst = [next(rr).rstrip() for x in range(4)]
    i1_lst = [next(fi).rstrip() for x in range(4)]
    i2_lst = [next(ri).rstrip() for x in range(4)]
    # find barcode id for I1 and I2 within barcode dict (=I1.id, I2.id) 
    bar1 = i1_lst[1]
    bar2 = revc(i2_lst[1])
    new_head_r1 = r1_lst[0] + " " + bar1 + "-" + bar2 
    new_head_r2 = r2_lst[0] + " " + bar1 + "-" + bar2 
    if bar1 not in bar_dict.values() or bar2 not in bar_dict.values():
        junk_fw.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
        junk_rv.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n")
    else:
        ## check for Q-scores
        q=0 # flag for bad qscore
        for basepos in range(len(i1_lst[3])):
            p=int(Bioinfo.convert_phred(i1_lst[3][basepos]))
            if p < minq:
                q = 1
        for basepos in range(len(i2_lst[3])):
            p=int(Bioinfo.convert_phred(i2_lst[3][basepos]))
            if p < minq:
                q = 1
        if q == 1: # if bad qscore, write to junk bins
            junk_fw.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
            junk_rv.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n")
        ## check for index swapping
        elif bar1 != bar2:
            swap_fw.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
            swap_rv.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n")            
        ## write everything else to sample bins
        else:
            index_id = [key for key, val in bar_dict.items() if val == bar1][0]
            f = [val for key, val in file_dict.items() if key == "FW_" + index_id][0]
            r = [val for key, val in file_dict.items() if key == "RV_" + index_id][0]
            f.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
            r.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n")
            
## close output files
swap_fw.close()
swap_rv.close() 
junk_fw.close()
junk_rv.close()
for file in file_dict.values():
    file.close()

## close input files
fr.close()
rr.close() 
fi.close()
ri.close()
