#!/usr/bin/python
# Jun-Hoe, Lee (2020)
# moidified from statsFasta.py to generate a basic seq length frequency histogram
# 1. get a distribution count of sequence length
# 2. generate a histogram of length distribution
# Usage: python  [input fasta file] [make output file]

import sys
import os
import argparse
from Bio.SeqUtils import GC
from Bio import SeqIO
import numpy as np

from matplotlib import pyplot as plt

######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# parser for input files and optional arguments
parser = errorDisplayParser(description='seqLengthDist: Get distribution of seq lengths in a \
                                         fasta file, then create a histogram')
parser.add_argument('inputSeqFile', action='store',
                     help='input fasta file, preferably in single line format')
parser.add_argument('outputHistFile', action='store',
                     help='output png file of histogram')
parser.add_argument('-bins', type=int, default=20,
                     help='number of bins, default=20. Usage: -bin 20')
parser.add_argument('-verbose', action='store_true', default=False,
                     help='turns on verbose mode. Usage: -verbose')  # verbose flag

args = parser.parse_args()

# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

####################################
def roundUp(x): # round up to next 100
    return x if x % 100 == 0 else x + 100 - x % 100


####################################
# A. read input file, get seq length distribution

try:
    seqs = [seq for seq in SeqIO.parse(args.inputSeqFile, "fasta")]
except:
    print("Incorrect input file")
    sys.exit(1)

lengths_list = [len(i.seq) for i in seqs]
if args.verbose:
    print("lengths_list :5 ", lengths_list[:5])

######################################
# B. generate histogram

# binsNum = np.linspace(math.ceil(min(lengths_list)),
#                    math.floor(max(lengths_list)),
#                    args.bins) # number of bins

plt.xlim([min(lengths_list)-5, max(lengths_list)+5])

plt.hist(lengths_list, bins=args.bins, alpha=0.5, edgecolor='black')
plt.xticks(np.arange(100 - roundUp(min(lengths_list)),
           roundUp(max(lengths_list)), 50))

plt.title('Distribution of sequence lengths')
plt.xlabel('Sequence length (bp)')
plt.ylabel('Number of sequences')

plt.savefig(args.outputHistFile)
