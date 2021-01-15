#!/usr/bin/python
# Jun-Hoe, Lee (2020)
# moidified from countGCperGenome.py to output more basic stats for a fasta file
# stats are numbe of seqs, average seq Length, longestSeq, shortestSeq, Gc%
# Usage: python countGCperGenome.py [input fasta file] [make output file]

import sys
import os
from Bio.SeqUtils import GC
from Bio import SeqIO
import numpy as np


##################################3
# function to turn on verbose mode
verbose = False
if verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

####################################


# A. get input file

try:
    fasta_file = sys.argv[1]
    verboseprint("fasta_file", fasta_file)
except:
    print("Incorrect input parameters")
    print("python countGCperGenome.py [input fasta file] [yes/no make output file]")
    sys.exit(1)


seqs = [seq for seq in SeqIO.parse(fasta_file, "fasta")]

lengths_list = [len(i.seq) for i in seqs]
# number of seqs
numSeqs = len(lengths_list)
# get average length
avgSeqLength = np.average(lengths_list)
# get longest and shortest seq
maxSeqLength = max(lengths_list)
minSeqLength = min(lengths_list)
## calculate GC%
gc_percent = np.average([GC(i.seq) for i in seqs])

# total sequnce length (not used for now)
total_size = np.sum(lengths_list)

print("fasta_file\tNumSeqs\tAvgSeqLength\tMaxSeqLength\tMinSeqLength\tGC_percent")
writeLine = "%s\t%d\t%0.2f\t%d\t%d\t%0.2f\n" % \
            (fasta_file, numSeqs, avgSeqLength, maxSeqLength, minSeqLength, gc_percent)
print(writeLine)

# check whether to create output file
if len(sys.argv) >  2:
    getInputFileName = os.path.basename(fasta_file)
    outputFileName = "out_" + getInputFileName
    verboseprint("make outputFileName", outputFileName)
    with open(outputFileName, 'w') as outputFile:
        outputFile.write(writeLine)
    outputFile.close()
