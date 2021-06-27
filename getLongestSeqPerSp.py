#/usr/bin/python3
# Jun-Hoe, Lee (2019)
# go through a fasta file of nucleotide or protein sequences to:
# 1. get longest seq
# input is a fasta file generated from e.g. EMBOSS that lists several ORF/prots per species
# ensure fasta file is already converted to 1 line

import os
import sys
import subprocess
import argparse
import configparser
import collections


from numpy import percentile
from Bio.Seq import Seq

######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# parser for input files and optional arguments

parser = errorDisplayParser(description='filterLengthRename: Go through a fasta file, rename the seqs, remove ambiguous residues, and/or filter each sequence by length and minimum number of species per gene ')
parser.add_argument('inputSeqFile', action='store',
                     help='input file of multiline fasta, in fasta format')
parser.add_argument('outputSeqFile', type=argparse.FileType('w'),
                     help='output fasta file of single line fasta')
parser.add_argument('separator', action="store", help='specify separator to get species name, \
                     e.g ">bo_Kamp_v1_2113_4_1_1", the separator is "_"')
parser.add_argument('uniqSpName_digit', action="store", type=int, help='specify the number of fragments to keep after \
                     split by separator, e.g. ">bo_Kamp_v1_2113_4_1_1", the number is 2 to keep "bo_Kamp" after split by "_"')
parser.add_argument('-verbose', action='store_true', default=False, help='turns on verbose mode. Usage: -verbose')  # verbose flag

args = parser.parse_args()


# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

###################################################
# function to get species name:
def getSpeciesName(split_symbol, sp_fragmentNum, givenLine):
    split_line = givenLine.split(splitter)
    if sp_fragmentNum == 1:
        speciesName = split_line[0]
    else:
        speciesName = split_symbol.join(map(str, split_line[:sp_fragmentNum]))
    return speciesName

######################################################
# A. read in files and find longest seq

splitter = args.separator
sp_digit = args.uniqSpName_digit

# create dicts
seqNameSeq_dict = collections.OrderedDict() # store the modified sequnece names and sequences in ordered dict
seqLength_list = [] # store sequence length
prevSp = "" # initiate empty string
uniqSeqName_list = [] # list to store uniq name of longest sequence

with open(args.inputSeqFile, "r") as readSeqFile:
    splitConvertStdOut = readSeqFile.readlines()

for line, lineNumber in zip(splitConvertStdOut, range(len(splitConvertStdOut))):
    line = line.rstrip('\n')
    if ">" in line: # is sequence names
        currentSp = getSpeciesName(args.separator, args.uniqSpName_digit, line)
        verboseprint("currentSp, line", currentSp, line)
        if currentSp != prevSp:
            seqNameSeq_dict[currentSp] = splitConvertStdOut[lineNumber + 1].rstrip('\n')  # first entry uses currentSp
            uniqSeqName_list.append(line)
            verboseprint("first new species added to dict ", line)
    else:
        seq_length = len(line) # length of sequence
        if seq_length > len(seqNameSeq_dict[currentSp]): # if the new sequence is longer
            seqNameSeq_dict[currentSp] = line # replaces old sequences
            uniqSeqName_list.pop()# remove last added
            prevSeqName = splitConvertStdOut[lineNumber - 1] # seq name from previous line
            uniqSeqName_list.append(prevSeqName.rstrip('\n'))
            verboseprint("len current seq lenght, prev seq length ", seq_length, len(seqNameSeq_dict[currentSp]))
    prevSp = currentSp
    verboseprint("prevSp, currentSp ", prevSp, currentSp)

verboseprint("uniqSeqName_list ", uniqSeqName_list)

if args.verbose:
    wrapNameSeq_dict = list(seqNameSeq_dict.items())
    print("seqNameSeq_dict", wrapNameSeq_dict[:5])


#####################################################
# B. write out to file

# write out the modified seq names seqs with filtering first:
newLine = ''

# the order in the dict should be the same as in uniqSeqName_list
seqCounter = 0
for seqHeader, longestSeq in seqNameSeq_dict.items():
    checkSp_name = getSpeciesName(args.separator, args.uniqSpName_digit, uniqSeqName_list[seqCounter])
    if seqHeader == checkSp_name:
        writeLine = uniqSeqName_list[seqCounter] + "\n" + longestSeq
#    writeLine = seqHeader + "\n" + longestSeq
        writeOutline = newLine + writeLine
        newLine = '\n'
        args.outputSeqFile.write(writeOutline)
        seqCounter += 1
    else:
        print("Error between unique species names and common species name of longest seq")
        sys.exit(2)
