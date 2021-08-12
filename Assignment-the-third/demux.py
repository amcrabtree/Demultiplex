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

# ## build this if you have time for error correction
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
## in bar_dict, value = index id, key = barcode
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
if f_reads.split(".")[-1] == "gz":
    ## open zipped input files
    fr = gzip.open(f_reads, "rt")
    rr = gzip.open(r_reads, "rt")
    fi = gzip.open(f_index, "rt")
    ri = gzip.open(r_index, "rt")
    ## count lines in file (will be same for all files)
    numrecs = int((max(i for i, _ in enumerate(ri))+1)/4)
    ri = gzip.open(r_index, "rt")
else:
    ## open unzipped input files
    fr = open(f_reads, "r")
    rr = open(r_reads, "r")
    fi = open(f_index, "r")
    ri = open(r_index, "r")
    ## count lines in file (will be same for all files)
    numrecs = int((max(i for i, _ in enumerate(ri))+1)/4)
    ri = open(r_index, "r")

## generate/open output files (sample bins [# of barcodes *2], 
## swapping bins [2], dust bins [2])
swap_fw = open (dir_out+"/FW_swap.fq", "a")
swap_rv = open (dir_out+"/RV_swap.fq", "a")
junk_fw = open (dir_out+"/FW_junk.fq", "a")
junk_rv = open (dir_out+"/RV_junk.fq", "a")

## make dictionary of output file info while opening them
## file_dict: key=FW/RV & index, value=file handle
file_dict = {}
for el in bar_dict:
    f = dir_out + "/FW_" + el + ".fq"
    r = dir_out + "/RV_" + el + ".fq"
    fwh = "FW_"+el
    rvh = "RV_"+el
    file_dict[fwh] = open(f, "w")
    file_dict[rvh] = open(r, "w")

## make dictionary of record counts to sum while binning
## count_dict: key=bin pair, value=number of record pairs
count_dict = {}
count_dict["junk"]=0
count_dict["swap"]=0
for el in bar_dict:
    count_dict[el]=0

## sort reads by indexes and write to output files
for i in range(numrecs):
    ## Loop for each fastq record in every read&index file at once
    r1_lst = [next(fr).rstrip() for x in range(4)]
    r2_lst = [next(rr).rstrip() for x in range(4)]
    i1_lst = [next(fi).rstrip() for x in range(4)]
    i2_lst = [next(ri).rstrip() for x in range(4)]
    bar1 = i1_lst[1]
    bar2 = revc(i2_lst[1])
    ## create new fastq headers for output files
    new_head_r1 = r1_lst[0] + " " + bar1 + "-" + bar2 
    new_head_r2 = r2_lst[0] + " " + bar1 + "-" + bar2 
    ## check if there are barcode sequences that are not in the barcode list (seq errors)
    if bar1 not in bar_dict.values() or bar2 not in bar_dict.values():
        junk_fw.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
        junk_rv.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n")
        count_dict["junk"] += 1
    else:
        ## check for Q-scores
        q=0 # flag for bad qscore
        for basepos in range(len(i1_lst[3])):  # check qscore of forward read index
            p=int(Bioinfo.convert_phred(i1_lst[3][basepos]))
            if p < minq:
                q = 1
        for basepos in range(len(i2_lst[3])):  # check qscore of reverse read index
            p=int(Bioinfo.convert_phred(i2_lst[3][basepos]))
            if p < minq:
                q = 1
        if q == 1:  # if bad qscore, write to junk bins
            junk_fw.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
            junk_rv.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n")
            count_dict["junk"] += 1
        ## check for index swapping
        elif bar1 != bar2:
            swap_fw.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
            swap_rv.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n") 
            count_dict["swap"] += 1           
        ## write everything else to sample bins
        else: 
            index_id = [key for key, val in bar_dict.items() if val == bar1][0]
            f = [val for key, val in file_dict.items() if key == "FW_" + index_id][0]
            r = [val for key, val in file_dict.items() if key == "RV_" + index_id][0]
            f.write(new_head_r1+"\n"+r1_lst[1]+"\n"+r1_lst[2]+"\n"+r1_lst[3]+"\n")
            r.write(new_head_r2+"\n"+r2_lst[1]+"\n"+r2_lst[2]+"\n"+r2_lst[3]+"\n")
            count_dict[index_id] += 1
            
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

## print count stats
for el in count_dict:
    print("Number of record pairs in", el, "bin:", count_dict[el])

########################## GRAPH ##############################

## print plot
import matplotlib.pyplot as plt
xval=[x for x in count_dict.keys()]
yval=[y for y in count_dict.values()]
plt.bar(xval, yval, color='blue')
plt.title("Number of Read Pairs Per Index")
plt.xlabel("Index/Bin")
plt.ylabel("Number of Read Pairs")
plt.savefig(dir_out+"/demux_hist.png", format="png")

########################## OUTPUT DATA ##############################

## write bin values to output file 
import sys
cmdline = (' '.join(sys.argv))
with open(dir_out+"/bin_counts.csv", "a") as mdrpt:
    mdrpt.write(cmdline + "\n")
    mdrpt.write("bin,rpcount\n")
    for el in count_dict:
        mdrpt.write(el + "," + str(count_dict[el]) + "\n")
